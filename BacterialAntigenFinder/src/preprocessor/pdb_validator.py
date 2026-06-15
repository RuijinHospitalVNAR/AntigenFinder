"""
PDB结构验证模块

验证PDB格式的蛋白质结构文件，提取结构信息。
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from dataclasses import dataclass, field

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue


@dataclass
class ProteinStructure:
    """蛋白质结构数据类"""
    id: str
    file_path: str
    chains: List[str]
    residue_count: int
    resolution: Optional[float] = None
    is_alphafold: bool = False
    sequence: str = ""
    plddt_scores: List[float] = field(default_factory=list)
    
    def get_chain_sequence(self, chain_id: str) -> str:
        """获取指定链的序列"""
        # 这里需要实际从结构中提取
        return self.sequence


class PdbValidator:
    """PDB文件验证器"""
    
    # 标准氨基酸三字母到单字母映射
    AA_3TO1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'SEC': 'U', 'PYL': 'O'  # 非标准氨基酸
    }
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.parser = PDBParser(QUIET=True)
    
    def validate(self, pdb_path: Union[str, Path]) -> Optional[ProteinStructure]:
        """
        验证单个PDB文件
        
        Args:
            pdb_path: PDB文件路径
        
        Returns:
            ProteinStructure对象，验证失败返回None
        """
        pdb_path = Path(pdb_path)
        
        if not pdb_path.exists():
            self.logger.error(f"PDB文件不存在: {pdb_path}")
            return None
        
        try:
            structure = self.parser.get_structure(pdb_path.stem, str(pdb_path))
            
            # 获取链信息
            chains = []
            residue_count = 0
            sequence = ""
            plddt_scores = []
            
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    chains.append(chain_id)
                    
                    # 提取序列和残基数
                    for residue in chain:
                        if residue.id[0] == ' ':  # 标准残基
                            residue_count += 1
                            resname = residue.resname
                            if resname in self.AA_3TO1:
                                sequence += self.AA_3TO1[resname]
                            else:
                                sequence += 'X'  # 未知残基
                            
                            # 提取pLDDT分数（如果是AlphaFold结构）
                            for atom in residue:
                                if atom.name == 'CA':
                                    plddt_scores.append(atom.bfactor)
                                    break
                break  # 只处理第一个模型
            
            # 判断是否为AlphaFold结构
            is_alphafold = self._is_alphafold_structure(pdb_path, plddt_scores)
            
            protein_struct = ProteinStructure(
                id=pdb_path.stem,
                file_path=str(pdb_path),
                chains=chains,
                residue_count=residue_count,
                is_alphafold=is_alphafold,
                sequence=sequence,
                plddt_scores=plddt_scores
            )
            
            # 验证结构完整性
            if not self._validate_structure(structure):
                self.logger.warning(f"结构 {pdb_path.stem} 可能不完整")
            
            return protein_struct
            
        except Exception as e:
            self.logger.error(f"解析PDB文件失败 {pdb_path}: {e}")
            return None
    
    def validate_directory(self, pdb_dir: Union[str, Path]) -> Dict[str, ProteinStructure]:
        """
        验证目录中的所有PDB文件
        
        Args:
            pdb_dir: PDB文件目录
        
        Returns:
            结构ID到ProteinStructure对象的字典
        """
        pdb_dir = Path(pdb_dir)
        
        if not pdb_dir.exists():
            raise FileNotFoundError(f"PDB目录不存在: {pdb_dir}")
        
        structures = {}
        pdb_files = list(pdb_dir.glob("*.pdb")) + list(pdb_dir.glob("*.PDB"))
        
        self.logger.info(f"在 {pdb_dir} 中找到 {len(pdb_files)} 个PDB文件")
        
        for pdb_file in pdb_files:
            struct = self.validate(pdb_file)
            if struct:
                structures[struct.id] = struct
        
        valid_count = len(structures)
        invalid_count = len(pdb_files) - valid_count
        
        if invalid_count > 0:
            self.logger.warning(f"{invalid_count} 个PDB文件验证失败")
        
        self.logger.info(f"成功验证 {valid_count} 个结构")
        
        return structures
    
    def _validate_structure(self, structure: Structure) -> bool:
        """
        验证结构完整性
        
        检查是否存在骨架原子 (N, CA, C)
        """
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':
                        continue
                    
                    # 检查骨架原子
                    atom_names = [atom.name for atom in residue]
                    if not all(atom in atom_names for atom in ['N', 'CA', 'C']):
                        return False
            break
        
        return True
    
    def _is_alphafold_structure(self, pdb_path: Path, 
                                  plddt_scores: List[float]) -> bool:
        """
        判断是否为AlphaFold预测结构
        
        基于文件名和B-factor范围判断
        """
        # 检查文件名
        name_lower = pdb_path.stem.lower()
        if 'alphafold' in name_lower or 'af-' in name_lower:
            return True
        
        # 检查B-factor范围（pLDDT在0-100之间）
        if plddt_scores:
            avg_plddt = sum(plddt_scores) / len(plddt_scores)
            # AlphaFold结构的pLDDT通常在0-100范围内
            if 0 < avg_plddt <= 100:
                # 进一步检查分布
                if all(0 <= score <= 100 for score in plddt_scores):
                    return True
        
        return False
    
    def extract_chain(self, pdb_path: Union[str, Path], 
                      chain_id: str,
                      output_path: Optional[Union[str, Path]] = None) -> str:
        """
        提取单条链到新的PDB文件
        
        Args:
            pdb_path: 原始PDB文件路径
            chain_id: 要提取的链ID
            output_path: 输出文件路径，默认为原文件名_链ID.pdb
        
        Returns:
            输出文件路径
        """
        pdb_path = Path(pdb_path)
        
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}_{chain_id}.pdb"
        else:
            output_path = Path(output_path)
        
        structure = self.parser.get_structure(pdb_path.stem, str(pdb_path))
        
        class ChainSelect(Select):
            def accept_chain(self, chain):
                return chain.id == chain_id
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), ChainSelect())
        
        return str(output_path)
    
    def get_structure_stats(self, structures: Dict[str, ProteinStructure]) -> dict:
        """
        获取结构统计信息
        
        Args:
            structures: 结构字典
        
        Returns:
            统计信息字典
        """
        if not structures:
            return {}
        
        residue_counts = [s.residue_count for s in structures.values()]
        alphafold_count = sum(1 for s in structures.values() if s.is_alphafold)
        
        return {
            'count': len(structures),
            'total_residues': sum(residue_counts),
            'min_residues': min(residue_counts),
            'max_residues': max(residue_counts),
            'avg_residues': sum(residue_counts) / len(residue_counts),
            'alphafold_count': alphafold_count,
            'experimental_count': len(structures) - alphafold_count
        }
