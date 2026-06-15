"""
DiscoTope-3.0 包装器

基于ESM-IF1的构象B细胞表位预测。
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


class DiscotopeWrapper(BasePredictor):
    """DiscoTope-3.0 预测器包装器"""
    
    def __init__(self, 
                 model_path: str,
                 env_name: str = 'discotope_env',
                 threshold: float = 0.90,
                 use_gpu: bool = True,
                 timeout: int = 3600,
                 structure_type: str = 'solved'):
        """
        初始化DiscoTope-3.0包装器
        
        Args:
            model_path: DiscoTope-3.0目录路径
            env_name: Conda环境名称
            threshold: 校准分数阈值（默认0.90）
            use_gpu: 是否使用GPU
            timeout: 超时时间
            structure_type: 结构类型 ('solved' 或 'alphafold')
        """
        super().__init__(model_path, env_name, threshold, use_gpu, timeout)
        self.structure_type = structure_type
        self.script_path = self.model_path / 'discotope3' / 'main.py'
        self.models_dir = self.model_path / 'models'
    
    @property
    def name(self) -> str:
        return 'discotope'
    
    @property
    def input_type(self) -> str:
        return 'structure'
    
    def predict(self, 
                sequences: Dict[str, ProteinSequence],
                structures: Optional[Dict[str, ProteinStructure]] = None,
                organism_type: str = 'gram-') -> PredictionResult:
        """
        运行DiscoTope-3.0预测
        
        Args:
            sequences: 序列字典（用于ID映射）
            structures: 结构字典
            organism_type: 不使用
        
        Returns:
            预测结果
        """
        if not structures:
            self.logger.error("DiscoTope-3.0需要结构输入")
            return PredictionResult(
                predictor_name=self.name,
                predictions=[],
                metadata={'error': 'No structures provided'}
            )
        
        # 创建临时目录
        temp_dir = tempfile.mkdtemp(prefix='discotope_')
        pdb_dir = os.path.join(temp_dir, 'pdbs')
        output_dir = os.path.join(temp_dir, 'output')
        os.makedirs(pdb_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # 复制PDB文件到临时目录
            for struct_id, struct in structures.items():
                src_path = struct.file_path
                dst_path = os.path.join(pdb_dir, f"{struct_id}.pdb")
                shutil.copy(src_path, dst_path)
            
            # 确定结构类型
            struc_type = self._determine_structure_type(structures)
            
            # 构建命令参数
            args = [
                '--pdb_dir', pdb_dir,
                '--out_dir', output_dir,
                '--struc_type', struc_type,
                '--calibrated_score_epi_threshold', str(self.threshold)
            ]
            
            if not self.use_gpu:
                args.append('--cpu_only')
            
            # 添加模型目录
            if self.models_dir.exists():
                args.extend(['--models_dir', str(self.models_dir)])
            
            self.logger.info(f"运行DiscoTope-3.0，输入 {len(structures)} 个结构")
            
            # 运行预测
            result = self.run_in_env(
                str(self.script_path),
                args,
                working_dir=str(self.model_path)
            )
            
            if result.returncode != 0:
                self.logger.error(f"DiscoTope-3.0运行失败: {result.stderr}")
                return PredictionResult(
                    predictor_name=self.name,
                    predictions=[],
                    metadata={'error': result.stderr}
                )
            
            # 解析输出
            return self._parse_output(output_dir)
            
        finally:
            # 清理临时目录
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                self.logger.warning(f"清理临时目录失败: {e}")
    
    def _determine_structure_type(self, 
                                    structures: Dict[str, ProteinStructure]) -> str:
        """
        确定结构类型
        
        如果大多数是AlphaFold结构，返回'alphafold'，否则返回'solved'
        """
        if not structures:
            return self.structure_type
        
        alphafold_count = sum(1 for s in structures.values() if s.is_alphafold)
        
        if alphafold_count > len(structures) / 2:
            return 'alphafold'
        
        return 'solved'
    
    def _parse_output(self, output_dir: str) -> PredictionResult:
        """
        解析DiscoTope-3.0输出
        
        Args:
            output_dir: 输出目录
        
        Returns:
            预测结果
        """
        predictions = []
        protein_scores = {}
        
        # 查找所有CSV输出文件
        output_path = Path(output_dir)
        
        # DiscoTope输出结构: output_dir/pdb_id/output/pdb_id_chain_discotope3.csv
        csv_files = list(output_path.rglob('*_discotope3.csv'))
        
        for csv_file in csv_files:
            try:
                df = pd.read_csv(csv_file)
                
                # 提取蛋白质ID
                protein_id = csv_file.stem.replace('_discotope3', '')
                
                scores_list = []
                
                for _, row in df.iterrows():
                    residue_id = int(row.get('res_id', 0))
                    residue_name = str(row.get('residue', 'X'))
                    
                    # DiscoTope-3.0分数
                    score = float(row.get('DiscoTope-3.0_score', 0))
                    
                    # 校准分数（如果有）
                    calibrated_score = float(row.get('calibrated_score', score))
                    
                    # 是否为表位
                    is_epitope = row.get('epitope', False)
                    if isinstance(is_epitope, str):
                        is_epitope = is_epitope.lower() == 'true'
                    
                    # RSA和pLDDT
                    rsa = float(row.get('rsa', 0))
                    plddt = float(row.get('pLDDTs', 100))
                    
                    predictions.append(EpitopePrediction(
                        protein_id=protein_id,
                        residue_id=residue_id,
                        residue_name=residue_name,
                        score=score,
                        is_epitope=is_epitope,
                        confidence=calibrated_score,
                        additional_info={
                            'calibrated_score': calibrated_score,
                            'rsa': rsa,
                            'plddt': plddt,
                            'predictor': 'discotope'
                        }
                    ))
                    
                    scores_list.append(score)
                
                # 计算蛋白质级别分数
                if scores_list:
                    protein_scores[protein_id] = sum(scores_list) / len(scores_list)
                    
            except Exception as e:
                self.logger.error(f"解析DiscoTope输出文件 {csv_file} 失败: {e}")
        
        self.logger.info(f"DiscoTope-3.0预测完成，{len(predictions)} 个残基")
        
        return PredictionResult(
            predictor_name=self.name,
            predictions=predictions,
            protein_scores=protein_scores,
            metadata={
                'threshold': self.threshold,
                'structure_type': self.structure_type,
                'total_proteins': len(csv_files),
                'total_residues': len(predictions)
            }
        )
