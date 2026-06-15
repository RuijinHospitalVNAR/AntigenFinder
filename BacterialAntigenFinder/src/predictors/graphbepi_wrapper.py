"""
GraphBepi 包装器

基于图神经网络的B细胞表位预测。
"""

import os
import tempfile
import shutil
from pathlib import Path
from typing import Dict, Optional
import pandas as pd

from .base_predictor import (
    BasePredictor, PredictionResult, EpitopePrediction
)
from ..preprocessor.fasta_parser import ProteinSequence
from ..preprocessor.pdb_validator import ProteinStructure


class GraphbepiWrapper(BasePredictor):
    """GraphBepi 预测器包装器"""
    
    def __init__(self, 
                 model_path: str,
                 env_name: str = 'graphbepi_env',
                 threshold: float = 0.1763,
                 use_gpu: bool = True,
                 timeout: int = 3600):
        """
        初始化GraphBepi包装器
        
        Args:
            model_path: GraphBepi目录路径
            env_name: Conda环境名称
            threshold: 预测阈值（默认0.1763）
            use_gpu: 是否使用GPU
            timeout: 超时时间
        """
        super().__init__(model_path, env_name, threshold, use_gpu, timeout)
        self.script_path = self.model_path / 'test.py'
    
    @property
    def name(self) -> str:
        return 'graphbepi'
    
    @property
    def input_type(self) -> str:
        return 'structure'  # 也支持序列（通过ESMFold预测结构）
    
    def predict(self, 
                sequences: Dict[str, ProteinSequence],
                structures: Optional[Dict[str, ProteinStructure]] = None,
                organism_type: str = 'gram-') -> PredictionResult:
        """
        运行GraphBepi预测
        
        Args:
            sequences: 序列字典
            structures: 结构字典
            organism_type: 不使用
        
        Returns:
            预测结果
        """
        # 创建临时目录
        temp_dir = tempfile.mkdtemp(prefix='graphbepi_')
        output_dir = os.path.join(temp_dir, 'output')
        os.makedirs(output_dir, exist_ok=True)
        
        all_predictions = []
        protein_scores = {}
        
        try:
            if structures:
                # 使用PDB文件进行预测
                for struct_id, struct in structures.items():
                    try:
                        preds = self._predict_single_pdb(
                            struct.file_path, 
                            output_dir,
                            struct_id
                        )
                        all_predictions.extend(preds)
                        
                        # 计算蛋白质分数
                        if preds:
                            scores = [p.score for p in preds]
                            protein_scores[struct_id] = sum(scores) / len(scores)
                            
                    except Exception as e:
                        self.logger.error(f"GraphBepi预测 {struct_id} 失败: {e}")
            else:
                # 使用FASTA序列（通过ESMFold预测结构）
                self.logger.warning("GraphBepi使用FASTA序列模式（ESMFold预测结构）")
                
                temp_fasta = os.path.join(temp_dir, 'input.fasta')
                with open(temp_fasta, 'w') as f:
                    for seq in sequences.values():
                        f.write(f">{seq.id}\n{seq.sequence}\n")
                
                preds = self._predict_from_fasta(temp_fasta, output_dir)
                all_predictions.extend(preds)
            
            self.logger.info(f"GraphBepi预测完成，{len(all_predictions)} 个残基")
            
            return PredictionResult(
                predictor_name=self.name,
                predictions=all_predictions,
                protein_scores=protein_scores,
                metadata={
                    'threshold': self.threshold,
                    'total_residues': len(all_predictions)
                }
            )
            
        finally:
            # 清理临时目录
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                self.logger.warning(f"清理临时目录失败: {e}")
    
    def _predict_single_pdb(self, pdb_path: str, 
                            output_dir: str,
                            protein_id: str) -> list:
        """
        对单个PDB文件进行预测
        """
        args = [
            '-i', pdb_path,
            '-p',  # PDB模式
            '-o', output_dir,
            '-t', str(self.threshold),
            '--gpu', '0' if self.use_gpu else '-1'
        ]
        
        result = self.run_in_env(
            str(self.script_path),
            args,
            working_dir=str(self.model_path)
        )
        
        if result.returncode != 0:
            self.logger.warning(f"GraphBepi预测失败: {result.stderr[:200]}")
            return []
        
        # 解析输出
        return self._parse_single_output(output_dir, protein_id)
    
    def _predict_from_fasta(self, fasta_path: str, output_dir: str) -> list:
        """
        从FASTA序列预测（使用ESMFold生成结构）
        """
        args = [
            '-i', fasta_path,
            '-f',  # FASTA模式
            '-o', output_dir,
            '-t', str(self.threshold),
            '--gpu', '0' if self.use_gpu else '-1'
        ]
        
        result = self.run_in_env(
            str(self.script_path),
            args,
            working_dir=str(self.model_path)
        )
        
        if result.returncode != 0:
            self.logger.warning(f"GraphBepi FASTA预测失败: {result.stderr[:200]}")
            return []
        
        # 解析所有输出文件
        predictions = []
        output_path = Path(output_dir)
        
        for csv_file in output_path.glob('*.csv'):
            protein_id = csv_file.stem
            predictions.extend(self._parse_single_output(output_dir, protein_id))
        
        return predictions
    
    def _parse_single_output(self, output_dir: str, protein_id: str) -> list:
        """
        解析单个预测输出
        """
        predictions = []
        
        # GraphBepi输出格式: output_dir/protein_id.csv
        csv_path = os.path.join(output_dir, f"{protein_id}.csv")
        
        if not os.path.exists(csv_path):
            # 尝试其他可能的文件名
            output_path = Path(output_dir)
            csv_files = list(output_path.glob('*.csv'))
            if csv_files:
                csv_path = str(csv_files[0])
            else:
                return predictions
        
        try:
            df = pd.read_csv(csv_path)
            
            for _, row in df.iterrows():
                # GraphBepi输出: resn, score, is epitope
                residue_info = str(row.get('resn', ''))
                score = float(row.get('score', 0))
                is_epitope = int(row.get('is epitope', 0)) == 1
                
                # 解析残基信息
                if residue_info:
                    parts = residue_info.split('_') if '_' in residue_info else [residue_info]
                    residue_name = parts[0] if parts else 'X'
                    try:
                        residue_id = int(parts[-1]) if len(parts) > 1 else 0
                    except ValueError:
                        residue_id = 0
                else:
                    residue_name = 'X'
                    residue_id = 0
                
                predictions.append(EpitopePrediction(
                    protein_id=protein_id,
                    residue_id=residue_id,
                    residue_name=residue_name,
                    score=score,
                    is_epitope=is_epitope,
                    confidence=score,
                    additional_info={'predictor': 'graphbepi'}
                ))
                
        except Exception as e:
            self.logger.error(f"解析GraphBepi输出失败 {csv_path}: {e}")
        
        return predictions
