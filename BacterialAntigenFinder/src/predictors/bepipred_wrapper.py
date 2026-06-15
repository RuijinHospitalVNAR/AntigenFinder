"""
BepiPred-3.0 包装器

基于ESM-2的B细胞表位预测。
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


class BepipredWrapper(BasePredictor):
    """BepiPred-3.0 预测器包装器"""
    
    def __init__(self, 
                 model_path: str,
                 env_name: str = 'bepipred_env',
                 threshold: float = 0.1512,
                 use_gpu: bool = True,
                 timeout: int = 3600,
                 prediction_mode: str = 'vt_pred',
                 top_percent: int = 20):
        """
        初始化BepiPred-3.0包装器
        
        Args:
            model_path: BepiPred-3.0目录路径
            env_name: Conda环境名称
            threshold: 预测阈值（默认0.1512）
            use_gpu: 是否使用GPU
            timeout: 超时时间
            prediction_mode: 预测模式 ('vt_pred' 或 'mjv_pred')
            top_percent: 输出前N%候选
        """
        super().__init__(model_path, env_name, threshold, use_gpu, timeout)
        self.prediction_mode = prediction_mode
        self.top_percent = top_percent
        self.script_path = self.model_path / 'bepipred3_CLI.py'
    
    @property
    def name(self) -> str:
        return 'bepipred'
    
    @property
    def input_type(self) -> str:
        return 'sequence'
    
    def predict(self, 
                sequences: Dict[str, ProteinSequence],
                structures: Optional[Dict[str, ProteinStructure]] = None,
                organism_type: str = 'gram-') -> PredictionResult:
        """
        运行BepiPred-3.0预测
        
        Args:
            sequences: 序列字典
            structures: 不使用
            organism_type: 不使用
        
        Returns:
            预测结果
        """
        if not self.validate_input(sequences):
            return PredictionResult(
                predictor_name=self.name,
                predictions=[],
                metadata={'error': 'Invalid input'}
            )
        
        # 创建临时目录
        temp_dir = tempfile.mkdtemp(prefix='bepipred_')
        temp_fasta = os.path.join(temp_dir, 'input.fasta')
        output_dir = os.path.join(temp_dir, 'output')
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # 写入临时FASTA文件
            with open(temp_fasta, 'w') as f:
                for seq in sequences.values():
                    f.write(f">{seq.id}\n{seq.sequence}\n")
            
            # 构建命令参数
            args = [
                '-i', temp_fasta,
                '-o', output_dir,
                '-pred', self.prediction_mode,
                '-t', str(self.threshold),
                '-top', str(self.top_percent)
            ]
            
            self.logger.info(f"运行BepiPred-3.0，输入 {len(sequences)} 个序列")
            
            # 运行预测
            result = self.run_in_env(
                str(self.script_path),
                args,
                working_dir=str(self.model_path)
            )
            
            if result.returncode != 0:
                self.logger.error(f"BepiPred-3.0运行失败: {result.stderr}")
                return PredictionResult(
                    predictor_name=self.name,
                    predictions=[],
                    metadata={'error': result.stderr}
                )
            
            # 解析输出
            return self._parse_output(output_dir, sequences)
            
        finally:
            # 清理临时目录
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                self.logger.warning(f"清理临时目录失败: {e}")
    
    def _parse_output(self, output_dir: str, 
                      sequences: Dict[str, ProteinSequence]) -> PredictionResult:
        """
        解析BepiPred-3.0输出
        
        Args:
            output_dir: 输出目录
            sequences: 原始序列字典
        
        Returns:
            预测结果
        """
        predictions = []
        protein_scores = {}
        
        # 解析raw_output.csv
        raw_output_path = os.path.join(output_dir, 'raw_output.csv')
        
        if os.path.exists(raw_output_path):
            try:
                df = pd.read_csv(raw_output_path)
                
                # BepiPred输出格式: protein_id, position, residue, score, ...
                for _, row in df.iterrows():
                    protein_id = str(row.get('ID', row.get('protein_id', '')))
                    residue_id = int(row.get('Position', row.get('position', 0)))
                    residue_name = str(row.get('Residue', row.get('residue', 'X')))
                    score = float(row.get('BepiPred-3.0 score', row.get('score', 0)))
                    
                    # 使用滚动平均分数作为线性表位分数
                    linear_score = float(row.get('BepiPred-3.0 linear epitope score', score))
                    
                    is_epitope = score >= self.threshold
                    
                    predictions.append(EpitopePrediction(
                        protein_id=protein_id,
                        residue_id=residue_id,
                        residue_name=residue_name,
                        score=score,
                        is_epitope=is_epitope,
                        confidence=score,
                        additional_info={
                            'linear_score': linear_score,
                            'predictor': 'bepipred'
                        }
                    ))
                    
                    # 累计蛋白质分数
                    if protein_id not in protein_scores:
                        protein_scores[protein_id] = []
                    protein_scores[protein_id].append(score)
                
            except Exception as e:
                self.logger.error(f"解析BepiPred输出失败: {e}")
        else:
            # 尝试解析FASTA输出
            fasta_output = os.path.join(output_dir, 'Bcell_epitope_preds.fasta')
            if os.path.exists(fasta_output):
                predictions = self._parse_fasta_output(fasta_output, sequences)
        
        # 计算蛋白质级别平均分数
        protein_avg_scores = {
            pid: sum(scores) / len(scores) 
            for pid, scores in protein_scores.items()
        }
        
        self.logger.info(f"BepiPred-3.0预测完成，{len(predictions)} 个残基")
        
        return PredictionResult(
            predictor_name=self.name,
            predictions=predictions,
            protein_scores=protein_avg_scores,
            metadata={
                'threshold': self.threshold,
                'prediction_mode': self.prediction_mode,
                'total_proteins': len(sequences),
                'total_residues': len(predictions)
            }
        )
    
    def _parse_fasta_output(self, fasta_path: str,
                            sequences: Dict[str, ProteinSequence]) -> list:
        """
        解析FASTA格式输出
        
        大写字母表示表位，小写字母表示非表位
        """
        predictions = []
        
        try:
            current_id = None
            current_seq = ""
            
            with open(fasta_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id and current_seq:
                            predictions.extend(
                                self._parse_epitope_sequence(current_id, current_seq)
                            )
                        current_id = line[1:].split()[0]
                        current_seq = ""
                    else:
                        current_seq += line
                
                # 处理最后一个序列
                if current_id and current_seq:
                    predictions.extend(
                        self._parse_epitope_sequence(current_id, current_seq)
                    )
                    
        except Exception as e:
            self.logger.error(f"解析FASTA输出失败: {e}")
        
        return predictions
    
    def _parse_epitope_sequence(self, protein_id: str, seq: str) -> list:
        """
        解析表位标记序列
        
        大写=表位，小写=非表位
        """
        predictions = []
        
        for i, char in enumerate(seq):
            is_epitope = char.isupper()
            residue_name = char.upper()
            
            # 简单的评分：表位为1.0，非表位为0.0
            score = 1.0 if is_epitope else 0.0
            
            predictions.append(EpitopePrediction(
                protein_id=protein_id,
                residue_id=i + 1,
                residue_name=residue_name,
                score=score,
                is_epitope=is_epitope,
                confidence=score,
                additional_info={'predictor': 'bepipred'}
            ))
        
        return predictions
