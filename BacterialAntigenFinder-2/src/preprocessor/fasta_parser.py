"""
FASTA序列解析模块

解析FASTA格式的蛋白质序列文件，提取序列ID和序列内容。
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Union
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Seq import Seq


@dataclass
class ProteinSequence:
    """蛋白质序列数据类"""
    id: str
    name: str
    description: str
    sequence: str
    length: int
    
    def __post_init__(self):
        self.length = len(self.sequence)
    
    def to_fasta(self) -> str:
        """转换为FASTA格式字符串"""
        return f">{self.id} {self.description}\n{self.sequence}\n"
    
    def get_subsequence(self, start: int, end: int) -> str:
        """获取子序列"""
        return self.sequence[start:end]
    
    def validate(self) -> bool:
        """验证序列是否为有效的蛋白质序列"""
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        return all(aa.upper() in valid_aa for aa in self.sequence if aa != '*')


class FastaParser:
    """FASTA文件解析器"""
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def parse(self, fasta_path: Union[str, Path]) -> Dict[str, ProteinSequence]:
        """
        解析FASTA文件
        
        Args:
            fasta_path: FASTA文件路径
        
        Returns:
            序列ID到ProteinSequence对象的字典
        """
        fasta_path = Path(fasta_path)
        
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA文件不存在: {fasta_path}")
        
        sequences = {}
        
        try:
            for record in SeqIO.parse(str(fasta_path), "fasta"):
                seq_id = self._normalize_id(record.id)
                
                # 提取名称和描述
                description = record.description
                name = description.split()[0] if description else seq_id
                
                # 创建序列对象
                protein_seq = ProteinSequence(
                    id=seq_id,
                    name=name,
                    description=description,
                    sequence=str(record.seq).upper(),
                    length=len(record.seq)
                )
                
                # 验证序列
                if not protein_seq.validate():
                    self.logger.warning(f"序列 {seq_id} 包含非标准氨基酸")
                
                sequences[seq_id] = protein_seq
                
            self.logger.info(f"成功解析 {len(sequences)} 个序列从 {fasta_path}")
            
        except Exception as e:
            self.logger.error(f"解析FASTA文件失败: {e}")
            raise
        
        return sequences
    
    def parse_string(self, fasta_string: str) -> Dict[str, ProteinSequence]:
        """
        解析FASTA格式字符串
        
        Args:
            fasta_string: FASTA格式的字符串
        
        Returns:
            序列ID到ProteinSequence对象的字典
        """
        from io import StringIO
        
        sequences = {}
        
        for record in SeqIO.parse(StringIO(fasta_string), "fasta"):
            seq_id = self._normalize_id(record.id)
            
            protein_seq = ProteinSequence(
                id=seq_id,
                name=record.id,
                description=record.description,
                sequence=str(record.seq).upper(),
                length=len(record.seq)
            )
            
            sequences[seq_id] = protein_seq
        
        return sequences
    
    def write(self, sequences: Dict[str, ProteinSequence], 
              output_path: Union[str, Path]) -> str:
        """
        将序列写入FASTA文件
        
        Args:
            sequences: 序列字典
            output_path: 输出文件路径
        
        Returns:
            输出文件路径
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            for seq in sequences.values():
                f.write(seq.to_fasta())
        
        self.logger.info(f"已写入 {len(sequences)} 个序列到 {output_path}")
        return str(output_path)
    
    def filter_by_length(self, sequences: Dict[str, ProteinSequence],
                         min_length: int = 50,
                         max_length: int = 5000) -> Dict[str, ProteinSequence]:
        """
        根据长度过滤序列
        
        Args:
            sequences: 序列字典
            min_length: 最小长度
            max_length: 最大长度
        
        Returns:
            过滤后的序列字典
        """
        filtered = {
            seq_id: seq 
            for seq_id, seq in sequences.items()
            if min_length <= seq.length <= max_length
        }
        
        filtered_count = len(sequences) - len(filtered)
        if filtered_count > 0:
            self.logger.info(f"根据长度过滤掉 {filtered_count} 个序列")
        
        return filtered
    
    def _normalize_id(self, seq_id: str) -> str:
        """
        规范化序列ID
        
        移除特殊字符，保留字母数字和下划线
        """
        # 保留基本ID部分
        normalized = re.sub(r'[^\w\-\.]', '_', seq_id)
        return normalized
    
    def get_sequence_stats(self, sequences: Dict[str, ProteinSequence]) -> dict:
        """
        获取序列统计信息
        
        Args:
            sequences: 序列字典
        
        Returns:
            统计信息字典
        """
        if not sequences:
            return {}
        
        lengths = [seq.length for seq in sequences.values()]
        
        return {
            'count': len(sequences),
            'total_residues': sum(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'avg_length': sum(lengths) / len(lengths),
            'median_length': sorted(lengths)[len(lengths) // 2]
        }
