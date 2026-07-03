"""
EpiGraph 包装器

基于图注意力网络的B细胞表位预测。
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


class EpigraphWrapper(BasePredictor):
    """EpiGraph 预测器包装器"""
    
    def __init__(self, 
                 model_path: str,
                 env_name: str = 'epigraph_env',
                 threshold: float = 0.1481,
                 use_gpu: bool = True,
                 timeout: int = 3600,
                 kfold: int = 10):
        """
        初始化EpiGraph包装器
        
        Args:
            model_path: EpiGraph目录路径
            env_name: Conda环境名称
            threshold: 分类阈值（默认0.1481）
            use_gpu: 是否使用GPU
            timeout: 超时时间
            kfold: 集成模型数量（默认10）
        """
        super().__init__(model_path, env_name, threshold, use_gpu, timeout)
        self.kfold = kfold
        self.inference_script = self.model_path / 'inference.py'
        self.custom_inference_script = self.model_path / 'inference_customPDB.py'
    
    @property
    def name(self) -> str:
        return 'epigraph'
    
    @property
    def input_type(self) -> str:
        return 'structure'
    
    def predict(self, 
                sequences: Dict[str, ProteinSequence],
                structures: Optional[Dict[str, ProteinStructure]] = None,
                organism_type: str = 'gram-') -> PredictionResult:
        """
        运行EpiGraph预测
        
        Args:
            sequences: 序列字典
            structures: 结构字典
            organism_type: 不使用
        
        Returns:
            预测结果
        """
        if not structures:
            self.logger.error("EpiGraph需要结构输入")
            return PredictionResult(
                predictor_name=self.name,
                predictions=[],
                metadata={'error': 'No structures provided'}
            )
        
        # 创建临时目录
        temp_dir = tempfile.mkdtemp(prefix='epigraph_')
        pdb_dir = os.path.join(temp_dir, 'Custom_PDB')
        result_dir = os.path.join(temp_dir, 'Result')
        os.makedirs(pdb_dir, exist_ok=True)
        os.makedirs(result_dir, exist_ok=True)
        
        all_predictions = []
        protein_scores = {}
        
        try:
            # 复制PDB文件到临时目录
            for struct_id, struct in structures.items():
                src_path = struct.file_path
                dst_path = os.path.join(pdb_dir, f"{struct_id}.pdb")
                shutil.copy(src_path, dst_path)
            
            # 对每个结构进行预测
            for struct_id in structures.keys():
                try:
                    preds = self._predict_single(struct_id, pdb_dir, result_dir)
                    all_predictions.extend(preds)
                    
                    # 计算蛋白质分数
                    if preds:
                        scores = [p.score for p in preds]
                        protein_scores[struct_id] = sum(scores) / len(scores)
                        
                except Exception as e:
                    self.logger.error(f"EpiGraph预测 {struct_id} 失败: {e}")
            
            self.logger.info(f"EpiGraph预测完成，{len(all_predictions)} 个残基")
            
            return PredictionResult(
                predictor_name=self.name,
                predictions=all_predictions,
                protein_scores=protein_scores,
                metadata={
                    'threshold': self.threshold,
                    'kfold': self.kfold,
                    'total_residues': len(all_predictions)
                }
            )
            
        finally:
            # 清理临时目录
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                self.logger.warning(f"清理临时目录失败: {e}")
    
    def _predict_single(self, pdb_id: str, pdb_dir: str, result_dir: str) -> list:
        """
        对单个PDB进行预测
        """
        args = [
            '--pdb', pdb_id,
            '--pdb_path', pdb_dir,
            '--out_path', result_dir,
            '--kfold', str(self.kfold),
            '--classification_threshold', str(self.threshold),
            '--device', 'cuda' if self.use_gpu else 'cpu'
        ]
        
        result = self.run_in_env(
            str(self.custom_inference_script),
            args,
            working_dir=str(self.model_path)
        )
        
        if result.returncode != 0:
            self.logger.warning(f"EpiGraph预测 {pdb_id} 失败: {result.stderr[:200]}")
            return []
        
        # 解析输出
        return self._parse_output(result_dir, pdb_id)
    
    def _parse_output(self, result_dir: str, pdb_id: str) -> list:
        """
        解析EpiGraph输出
        
        输出格式: PDB, Residue, Score, Epitope, RSA
        """
        predictions = []
        
        # EpiGraph输出: Result/pdb_id.csv
        csv_path = os.path.join(result_dir, f"{pdb_id}.csv")
        
        if not os.path.exists(csv_path):
            # 尝试查找其他可能的输出文件
            result_path = Path(result_dir)
            csv_files = list(result_path.glob(f"*{pdb_id}*.csv"))
            if csv_files:
                csv_path = str(csv_files[0])
            else:
                self.logger.warning(f"找不到EpiGraph输出文件: {pdb_id}")
                return predictions
        
        try:
            df = pd.read_csv(csv_path)
            
            for _, row in df.iterrows():
                # 解析残基信息 (格式: chain:res_name:res_id)
                residue_str = str(row.get('Residue', ''))
                parts = residue_str.split(':')
                
                if len(parts) >= 3:
                    chain = parts[0]
                    residue_name = parts[1]
                    try:
                        residue_id = int(parts[2])
                    except ValueError:
                        residue_id = 0
                else:
                    chain = ''
                    residue_name = residue_str
                    residue_id = 0
                
                score = float(row.get('Score', 0))
                is_epitope = int(row.get('Epitope', 0)) == 1
                rsa = float(row.get('RSA', 0))
                
                # Use pdb_id directly so predictions align with other structure-based
                # predictors (e.g. GraphBepi) that use the structure filename as id.
                protein_id = pdb_id
                
                predictions.append(EpitopePrediction(
                    protein_id=protein_id,
                    residue_id=residue_id,
                    residue_name=residue_name,
                    score=score,
                    is_epitope=is_epitope,
                    confidence=score,
                    additional_info={
                        'chain': chain,
                        'rsa': rsa,
                        'predictor': 'epigraph'
                    }
                ))
                
        except Exception as e:
            self.logger.error(f"解析EpiGraph输出失败 {csv_path}: {e}")
        
        return predictions
