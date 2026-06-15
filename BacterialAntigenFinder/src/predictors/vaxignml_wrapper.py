"""
Vaxign-ML 包装器

保护性抗原预测和免疫原性评估（通过Docker调用）。
"""

import os
import tempfile
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional
import pandas as pd

from .base_predictor import (
    BasePredictor, PredictionResult, EpitopePrediction
)
from ..preprocessor.fasta_parser import ProteinSequence
from ..preprocessor.pdb_validator import ProteinStructure


class VaxignmlWrapper(BasePredictor):
    """Vaxign-ML 预测器包装器（通过Docker调用）"""
    
    DOCKER_IMAGE = 'e4ong1031/vaxign-ml:latest'
    
    def __init__(self, 
                 model_path: str,
                 env_name: str = '',  # Docker不需要conda环境
                 threshold: float = 0.5,
                 use_gpu: bool = False,  # Vaxign-ML不使用GPU
                 timeout: int = 3600,
                 use_docker: bool = True):
        """
        初始化Vaxign-ML包装器
        
        Args:
            model_path: Vaxign-ML目录路径
            env_name: 不使用
            threshold: 预测阈值
            use_gpu: 不使用
            timeout: 超时时间
            use_docker: 是否使用Docker（True）或本地运行
        """
        super().__init__(model_path, env_name, threshold, use_gpu, timeout)
        self.use_docker = use_docker
        self.script_path = self.model_path / 'VaxignML.py'
    
    @property
    def name(self) -> str:
        return 'vaxignml'
    
    @property
    def input_type(self) -> str:
        return 'sequence'
    
    def predict(self, 
                sequences: Dict[str, ProteinSequence],
                structures: Optional[Dict[str, ProteinStructure]] = None,
                organism_type: str = 'gram-') -> PredictionResult:
        """
        运行Vaxign-ML预测
        
        Args:
            sequences: 序列字典
            structures: 不使用
            organism_type: 生物类型 ('gram+', 'gram-', 'virus')
        
        Returns:
            预测结果
        """
        if not self.validate_input(sequences):
            return PredictionResult(
                predictor_name=self.name,
                predictions=[],
                metadata={'error': 'Invalid input'}
            )
        
        # 规范化生物类型
        org_type = self._normalize_organism_type(organism_type)
        
        # 创建临时目录
        temp_dir = tempfile.mkdtemp(prefix='vaxignml_')
        temp_fasta = os.path.join(temp_dir, 'input.fasta')
        output_dir = os.path.join(temp_dir, 'output')
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # 写入FASTA文件
            with open(temp_fasta, 'w') as f:
                for seq in sequences.values():
                    f.write(f">{seq.id}\n{seq.sequence}\n")
            
            self.logger.info(f"运行Vaxign-ML，输入 {len(sequences)} 个序列，类型 {org_type}")
            
            if self.use_docker:
                success = self._run_docker(temp_fasta, output_dir, org_type)
            else:
                success = self._run_local(temp_fasta, output_dir, org_type)
            
            if not success:
                return PredictionResult(
                    predictor_name=self.name,
                    predictions=[],
                    metadata={'error': 'Vaxign-ML execution failed'}
                )
            
            # 解析输出
            return self._parse_output(output_dir, sequences)
            
        finally:
            # 清理临时目录
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                self.logger.warning(f"清理临时目录失败: {e}")
    
    def _normalize_organism_type(self, organism_type: str) -> str:
        """规范化生物类型"""
        org_map = {
            'gram+': 'gram+',
            'gram-': 'gram-',
            'g+': 'gram+',
            'g-': 'gram-',
            'virus': 'virus',
            'v': 'virus'
        }
        return org_map.get(organism_type.lower(), 'gram-')
    
    def _run_docker(self, fasta_path: str, output_dir: str, org_type: str) -> bool:
        """通过Docker运行Vaxign-ML"""
        # 创建PSORTB结果目录
        psortb_dir = os.path.join(output_dir, '_FEATURE', 'PSORTB')
        os.makedirs(psortb_dir, exist_ok=True)
        
        cmd = [
            'docker', 'run', '--rm',
            '-v', f'{fasta_path}:{fasta_path}',
            '-v', f'{output_dir}:{output_dir}',
            '-v', f'{psortb_dir}:/tmp/results',
            self.DOCKER_IMAGE,
            'python3.6', 'VaxignML.py',
            '-i', fasta_path,
            '-o', output_dir,
            '-t', org_type
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            if result.returncode != 0:
                self.logger.error(f"Vaxign-ML Docker运行失败: {result.stderr}")
                return False
            
            return True
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"Vaxign-ML执行超时 ({self.timeout}秒)")
            return False
        except FileNotFoundError:
            self.logger.error("Docker未安装或不在PATH中")
            return False
        except Exception as e:
            self.logger.error(f"Vaxign-ML执行失败: {e}")
            return False
    
    def _run_local(self, fasta_path: str, output_dir: str, org_type: str) -> bool:
        """本地运行Vaxign-ML（需要依赖环境）"""
        args = [
            '-i', fasta_path,
            '-o', output_dir,
            '-t', org_type
        ]
        
        try:
            result = subprocess.run(
                ['python', str(self.script_path)] + args,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=str(self.model_path)
            )
            
            if result.returncode != 0:
                self.logger.error(f"Vaxign-ML本地运行失败: {result.stderr}")
                return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"Vaxign-ML本地执行失败: {e}")
            return False
    
    def _parse_output(self, output_dir: str, 
                      sequences: Dict[str, ProteinSequence]) -> PredictionResult:
        """
        解析Vaxign-ML输出
        
        输出格式: sample, prediction, protegenicity
        """
        predictions = []
        protein_scores = {}
        
        # 查找结果文件
        output_path = Path(output_dir)
        result_files = list(output_path.glob('*.result.tsv'))
        
        if not result_files:
            # 尝试其他可能的输出格式
            result_files = list(output_path.glob('*result*.csv')) + \
                          list(output_path.glob('*result*.tsv'))
        
        for result_file in result_files:
            try:
                # 尝试读取TSV或CSV
                if result_file.suffix == '.tsv':
                    df = pd.read_csv(result_file, sep='\t')
                else:
                    df = pd.read_csv(result_file)
                
                for _, row in df.iterrows():
                    protein_id = str(row.get('sample', row.get('protein_id', '')))
                    prediction = float(row.get('prediction', 0))
                    protegenicity = float(row.get('protegenicity', prediction * 100))
                    
                    # Vaxign-ML是蛋白质级别的预测，不是残基级别
                    is_protective = prediction >= self.threshold
                    
                    # 存储蛋白质级别分数
                    protein_scores[protein_id] = protegenicity
                    
                    # 为每个残基创建预测（使用蛋白质分数）
                    if protein_id in sequences:
                        seq = sequences[protein_id]
                        for i, aa in enumerate(seq.sequence):
                            predictions.append(EpitopePrediction(
                                protein_id=protein_id,
                                residue_id=i + 1,
                                residue_name=aa,
                                score=prediction,
                                is_epitope=is_protective,
                                confidence=protegenicity / 100,
                                additional_info={
                                    'protegenicity': protegenicity,
                                    'is_protective_antigen': is_protective,
                                    'predictor': 'vaxignml'
                                }
                            ))
                    
            except Exception as e:
                self.logger.error(f"解析Vaxign-ML输出失败 {result_file}: {e}")
        
        self.logger.info(f"Vaxign-ML预测完成，{len(protein_scores)} 个蛋白质")
        
        return PredictionResult(
            predictor_name=self.name,
            predictions=predictions,
            protein_scores=protein_scores,
            metadata={
                'threshold': self.threshold,
                'total_proteins': len(protein_scores),
                'protective_antigens': sum(1 for s in protein_scores.values() if s >= self.threshold * 100)
            }
        )
    
    def check_docker(self) -> bool:
        """检查Docker是否可用"""
        try:
            result = subprocess.run(
                ['docker', '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def pull_image(self) -> bool:
        """拉取Docker镜像"""
        try:
            self.logger.info(f"拉取Docker镜像: {self.DOCKER_IMAGE}")
            result = subprocess.run(
                ['docker', 'pull', self.DOCKER_IMAGE],
                capture_output=True,
                text=True,
                timeout=600
            )
            return result.returncode == 0
        except Exception as e:
            self.logger.error(f"拉取Docker镜像失败: {e}")
            return False
