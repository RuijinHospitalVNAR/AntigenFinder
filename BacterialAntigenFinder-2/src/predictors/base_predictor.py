"""
预测器基类模块

定义所有预测器的通用接口和环境调用机制。
"""

import logging
import os
import subprocess
import tempfile
import json
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
from dataclasses import dataclass, field
import pandas as pd

from ..preprocessor.fasta_parser import ProteinSequence
from ..preprocessor.pdb_validator import ProteinStructure


@dataclass
class EpitopePrediction:
    """表位预测结果数据类"""
    protein_id: str
    residue_id: int
    residue_name: str
    score: float
    is_epitope: bool
    confidence: float = 0.0
    additional_info: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PredictionResult:
    """预测结果汇总数据类"""
    predictor_name: str
    predictions: List[EpitopePrediction]
    protein_scores: Dict[str, float] = field(default_factory=dict)  # 蛋白质级别分数
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dataframe(self) -> pd.DataFrame:
        """转换为DataFrame格式"""
        data = []
        for pred in self.predictions:
            row = {
                'protein_id': pred.protein_id,
                'residue_id': pred.residue_id,
                'residue_name': pred.residue_name,
                'score': pred.score,
                'is_epitope': pred.is_epitope,
                'confidence': pred.confidence
            }
            row.update(pred.additional_info)
            data.append(row)
        
        return pd.DataFrame(data)
    
    def get_epitopes(self, protein_id: Optional[str] = None) -> List[EpitopePrediction]:
        """获取被预测为表位的残基"""
        epitopes = [p for p in self.predictions if p.is_epitope]
        if protein_id:
            epitopes = [p for p in epitopes if p.protein_id == protein_id]
        return epitopes
    
    def get_epitope_regions(self, protein_id: str, 
                            min_length: int = 5) -> List[tuple]:
        """
        获取连续的表位区域
        
        Args:
            protein_id: 蛋白质ID
            min_length: 最小区域长度
        
        Returns:
            [(start, end, avg_score), ...]
        """
        epitopes = self.get_epitopes(protein_id)
        if not epitopes:
            return []
        
        # 按残基ID排序
        epitopes.sort(key=lambda x: x.residue_id)
        
        regions = []
        current_start = epitopes[0].residue_id
        current_end = current_start
        scores = [epitopes[0].score]
        
        for epitope in epitopes[1:]:
            if epitope.residue_id == current_end + 1:
                current_end = epitope.residue_id
                scores.append(epitope.score)
            else:
                if current_end - current_start + 1 >= min_length:
                    avg_score = sum(scores) / len(scores)
                    regions.append((current_start, current_end, avg_score))
                current_start = epitope.residue_id
                current_end = current_start
                scores = [epitope.score]
        
        # 处理最后一个区域
        if current_end - current_start + 1 >= min_length:
            avg_score = sum(scores) / len(scores)
            regions.append((current_start, current_end, avg_score))
        
        return regions


class BasePredictor(ABC):
    """预测器基类"""
    
    def __init__(self, 
                 model_path: str,
                 env_name: str,
                 threshold: float = 0.5,
                 use_gpu: bool = True,
                 timeout: int = 3600):
        """
        初始化预测器
        
        Args:
            model_path: 模型目录路径
            env_name: Conda环境名称
            threshold: 预测阈值
            use_gpu: 是否使用GPU
            timeout: 超时时间（秒）
        """
        self.model_path = Path(model_path).resolve()
        self.env_name = env_name
        self.threshold = threshold
        self.use_gpu = use_gpu
        self.timeout = timeout
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # 验证模型路径
        if not self.model_path.exists():
            self.logger.warning(f"模型路径不存在: {self.model_path}")
    
    @property
    @abstractmethod
    def name(self) -> str:
        """预测器名称"""
        pass
    
    @property
    @abstractmethod
    def input_type(self) -> str:
        """输入类型: 'sequence', 'structure', 'both'"""
        pass
    
    @abstractmethod
    def predict(self, 
                sequences: Dict[str, ProteinSequence],
                structures: Optional[Dict[str, ProteinStructure]] = None,
                organism_type: str = 'gram-') -> PredictionResult:
        """
        运行预测
        
        Args:
            sequences: 序列字典
            structures: 结构字典（可选）
            organism_type: 生物类型
        
        Returns:
            预测结果
        """
        pass
    
    def run_in_env(self, script: str, args: List[str], 
                    working_dir: Optional[str] = None,
                    env_vars: Optional[Dict[str, str]] = None) -> subprocess.CompletedProcess:
        """
        在指定Conda环境中运行脚本
        
        Args:
            script: 脚本路径或命令
            args: 命令行参数列表
            working_dir: 工作目录
            env_vars: 额外的环境变量
        
        Returns:
            subprocess.CompletedProcess对象
        """
        # 构建命令
        cmd = ['conda', 'run', '-n', self.env_name, 'python', script] + args
        
        # 设置环境变量
        env = os.environ.copy()
        if env_vars:
            env.update(env_vars)
        
        # 设置CUDA可见性
        if not self.use_gpu:
            env['CUDA_VISIBLE_DEVICES'] = ''
        
        self.logger.debug(f"执行命令: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=working_dir or str(self.model_path),
                env=env
            )
            
            if result.returncode != 0:
                self.logger.error(f"命令执行失败: {result.stderr}")
            
            return result
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时 ({self.timeout}秒)")
            raise
        except Exception as e:
            self.logger.error(f"命令执行异常: {e}")
            raise
    
    def run_command(self, cmd: List[str], 
                     working_dir: Optional[str] = None,
                     use_conda: bool = True) -> subprocess.CompletedProcess:
        """
        运行命令
        
        Args:
            cmd: 命令列表
            working_dir: 工作目录
            use_conda: 是否使用conda环境
        
        Returns:
            subprocess.CompletedProcess对象
        """
        if use_conda:
            full_cmd = ['conda', 'run', '-n', self.env_name] + cmd
        else:
            full_cmd = cmd
        
        env = os.environ.copy()
        if not self.use_gpu:
            env['CUDA_VISIBLE_DEVICES'] = ''
        
        self.logger.debug(f"执行命令: {' '.join(full_cmd)}")
        
        return subprocess.run(
            full_cmd,
            capture_output=True,
            text=True,
            timeout=self.timeout,
            cwd=working_dir or str(self.model_path),
            env=env
        )
    
    def prepare_temp_fasta(self, sequences: Dict[str, ProteinSequence]) -> str:
        """
        准备临时FASTA文件
        
        Args:
            sequences: 序列字典
        
        Returns:
            临时文件路径
        """
        fd, temp_path = tempfile.mkstemp(suffix='.fasta')
        
        with os.fdopen(fd, 'w') as f:
            for seq in sequences.values():
                f.write(f">{seq.id}\n{seq.sequence}\n")
        
        return temp_path
    
    def cleanup_temp_files(self, *paths: str):
        """清理临时文件"""
        for path in paths:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except Exception as e:
                self.logger.warning(f"清理临时文件失败 {path}: {e}")
    
    def parse_output(self, output_path: str) -> PredictionResult:
        """
        解析预测输出（子类可覆盖）
        
        Args:
            output_path: 输出文件路径
        
        Returns:
            预测结果
        """
        raise NotImplementedError("子类需要实现parse_output方法")
    
    def validate_input(self, sequences: Dict[str, ProteinSequence],
                       structures: Optional[Dict[str, ProteinStructure]] = None) -> bool:
        """
        验证输入数据
        
        Args:
            sequences: 序列字典
            structures: 结构字典
        
        Returns:
            是否有效
        """
        if not sequences:
            self.logger.error("没有输入序列")
            return False
        
        if self.input_type == 'structure' and not structures:
            self.logger.error(f"{self.name} 需要结构输入")
            return False
        
        return True
    
    def check_environment(self) -> bool:
        """
        检查Conda环境是否可用
        
        Returns:
            环境是否可用
        """
        try:
            result = subprocess.run(
                ['conda', 'run', '-n', self.env_name, 'python', '--version'],
                capture_output=True,
                text=True,
                timeout=30
            )
            return result.returncode == 0
        except Exception as e:
            self.logger.error(f"检查环境失败: {e}")
            return False
