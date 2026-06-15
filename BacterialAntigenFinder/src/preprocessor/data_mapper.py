"""
序列-结构映射模块

将FASTA序列与PDB结构进行匹配和映射。
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass

from .fasta_parser import ProteinSequence
from .pdb_validator import ProteinStructure


@dataclass
class MappedData:
    """序列-结构映射数据类"""
    id: str
    sequence: ProteinSequence
    structure: Optional[ProteinStructure]
    mapping_quality: float  # 0-1，表示序列与结构的匹配程度
    aligned_regions: List[Tuple[int, int]]  # 匹配的区域列表
    
    @property
    def has_structure(self) -> bool:
        return self.structure is not None
    
    @property
    def sequence_length(self) -> int:
        return self.sequence.length
    
    @property
    def structure_residues(self) -> int:
        return self.structure.residue_count if self.structure else 0


class DataMapper:
    """序列-结构映射器"""
    
    def __init__(self, strict_mode: bool = False):
        """
        初始化映射器
        
        Args:
            strict_mode: 严格模式下，只返回完全匹配的映射
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.strict_mode = strict_mode
    
    def map(self, sequences: Dict[str, ProteinSequence],
            structures: Dict[str, ProteinStructure]) -> Dict[str, MappedData]:
        """
        映射序列和结构
        
        Args:
            sequences: 序列字典
            structures: 结构字典
        
        Returns:
            映射数据字典
        """
        mapped_data = {}
        
        self.logger.info(f"开始映射 {len(sequences)} 个序列和 {len(structures)} 个结构")
        
        for seq_id, sequence in sequences.items():
            # 尝试查找匹配的结构
            structure = self._find_matching_structure(seq_id, sequence, structures)
            
            if structure:
                # 计算映射质量
                quality, aligned_regions = self._calculate_mapping_quality(
                    sequence, structure
                )
                
                if self.strict_mode and quality < 0.95:
                    self.logger.warning(
                        f"序列 {seq_id} 与结构匹配质量低 ({quality:.2f})，跳过"
                    )
                    structure = None
                    quality = 0.0
                    aligned_regions = []
            else:
                quality = 0.0
                aligned_regions = []
            
            mapped_data[seq_id] = MappedData(
                id=seq_id,
                sequence=sequence,
                structure=structure,
                mapping_quality=quality,
                aligned_regions=aligned_regions
            )
        
        # 统计映射结果
        with_structure = sum(1 for m in mapped_data.values() if m.has_structure)
        without_structure = len(mapped_data) - with_structure
        
        self.logger.info(f"映射完成: {with_structure} 个有结构, {without_structure} 个无结构")
        
        return mapped_data
    
    def _find_matching_structure(self, seq_id: str, 
                                   sequence: ProteinSequence,
                                   structures: Dict[str, ProteinStructure]) -> Optional[ProteinStructure]:
        """
        查找匹配的结构
        
        尝试多种匹配策略：
        1. 精确ID匹配
        2. 模糊ID匹配
        3. 序列匹配
        """
        # 1. 精确ID匹配
        if seq_id in structures:
            return structures[seq_id]
        
        # 2. 规范化ID后匹配
        normalized_seq_id = self._normalize_id(seq_id)
        for struct_id, structure in structures.items():
            if self._normalize_id(struct_id) == normalized_seq_id:
                return structure
        
        # 3. 部分ID匹配
        for struct_id, structure in structures.items():
            # 检查是否包含
            if normalized_seq_id in self._normalize_id(struct_id):
                return structure
            if self._normalize_id(struct_id) in normalized_seq_id:
                return structure
        
        # 4. 序列匹配（如果结构中有序列信息）
        for struct_id, structure in structures.items():
            if structure.sequence:
                # 检查序列相似度
                similarity = self._calculate_sequence_similarity(
                    sequence.sequence, structure.sequence
                )
                if similarity > 0.95:
                    self.logger.info(f"通过序列匹配: {seq_id} -> {struct_id}")
                    return structure
        
        return None
    
    def _normalize_id(self, id_str: str) -> str:
        """规范化ID用于匹配"""
        # 转小写，移除常见后缀
        normalized = id_str.lower()
        
        # 移除链标识符
        if '_' in normalized and len(normalized.split('_')[-1]) == 1:
            normalized = '_'.join(normalized.split('_')[:-1])
        
        # 移除常见后缀
        for suffix in ['.pdb', '_model', '_alphafold', '_af']:
            if normalized.endswith(suffix):
                normalized = normalized[:-len(suffix)]
        
        return normalized
    
    def _calculate_mapping_quality(self, sequence: ProteinSequence,
                                     structure: ProteinStructure) -> Tuple[float, List[Tuple[int, int]]]:
        """
        计算序列与结构的映射质量
        
        Returns:
            (质量分数, 对齐区域列表)
        """
        if not structure.sequence:
            # 无法验证，假设完全匹配
            return 1.0, [(0, sequence.length)]
        
        # 计算序列相似度
        similarity = self._calculate_sequence_similarity(
            sequence.sequence, structure.sequence
        )
        
        # 查找对齐区域
        aligned_regions = self._find_aligned_regions(
            sequence.sequence, structure.sequence
        )
        
        return similarity, aligned_regions
    
    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """
        计算两个序列的相似度
        
        使用简单的匹配率计算
        """
        if not seq1 or not seq2:
            return 0.0
        
        # 使用较短序列的长度
        min_len = min(len(seq1), len(seq2))
        
        if min_len == 0:
            return 0.0
        
        # 计算匹配残基数
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        
        return matches / min_len
    
    def _find_aligned_regions(self, seq1: str, seq2: str) -> List[Tuple[int, int]]:
        """
        查找对齐区域
        
        返回(start, end)元组列表
        """
        regions = []
        
        min_len = min(len(seq1), len(seq2))
        
        in_region = False
        start = 0
        
        for i in range(min_len):
            if seq1[i] == seq2[i]:
                if not in_region:
                    start = i
                    in_region = True
            else:
                if in_region:
                    regions.append((start, i))
                    in_region = False
        
        if in_region:
            regions.append((start, min_len))
        
        return regions
    
    def get_sequences_without_structures(self, 
                                          mapped_data: Dict[str, MappedData]) -> List[str]:
        """
        获取没有结构的序列ID列表
        """
        return [seq_id for seq_id, data in mapped_data.items() if not data.has_structure]
    
    def get_mapping_stats(self, mapped_data: Dict[str, MappedData]) -> dict:
        """
        获取映射统计信息
        """
        with_structure = sum(1 for m in mapped_data.values() if m.has_structure)
        
        qualities = [m.mapping_quality for m in mapped_data.values() if m.has_structure]
        avg_quality = sum(qualities) / len(qualities) if qualities else 0.0
        
        return {
            'total': len(mapped_data),
            'with_structure': with_structure,
            'without_structure': len(mapped_data) - with_structure,
            'coverage': with_structure / len(mapped_data) if mapped_data else 0.0,
            'avg_mapping_quality': avg_quality
        }
