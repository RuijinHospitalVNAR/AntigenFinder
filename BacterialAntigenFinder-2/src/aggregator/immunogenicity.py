"""
免疫原性评估模块

使用VaxiJen 2.0算法进行蛋白质级别的免疫原性评估。

VaxiJen 2.0基于氨基酸z-scale理化性质和自相关系数(ACC)方法预测抗原性，
无需外部依赖，纯Python+numpy实现，可完全在Docker中运行。
"""

import logging
from typing import Dict, List, Optional
from dataclasses import dataclass, field
import pandas as pd

from ..preprocessor.fasta_parser import ProteinSequence
from .vaxijen_calculator import VaxiJenCalculator, VaxiJenResult


@dataclass
class ImmunogenicityResult:
    """免疫原性评估结果"""
    protein_id: str
    protegenicity_score: float  # 保护性分数（基于VaxiJen抗原性映射到0-100）
    is_protective_antigen: bool
    subcellular_location: str  # 亚细胞定位
    has_signal_peptide: bool
    transmembrane_regions: int
    antigenicity_score: float  # VaxiJen 2.0抗原性分数 (0-1)
    immunogenicity_rank: int  # 免疫原性排名
    recommendation: str  # 推荐等级: HIGH, MEDIUM, LOW
    features: Dict[str, float] = field(default_factory=dict)


class ImmunogenicityEvaluator:
    """
    免疫原性评估器

    使用VaxiJen 2.0算法评估蛋白质抗原性，
    将抗原性分数映射为保护性分数用于候选排序。
    """
    
    # 优选亚细胞定位
    PREFERRED_LOCATIONS = [
        'OuterMembrane',
        'Extracellular',
        'Periplasmic',
        'CellWall',
        'Secreted',
        'Fimbrial',
        'Surface'
    ]
    
    def __init__(self, 
                 protegenicity_threshold: float = 50.0,
                 antigenicity_threshold: float = 0.4,
                 organism_type: str = 'bacteria',
                 use_vaxijen: bool = True):
        """
        初始化评估器
        
        Args:
            protegenicity_threshold: 保护性抗原阈值（百分位）
            antigenicity_threshold: 抗原性阈值
            organism_type: 生物类型 (用于VaxiJen模型选择)
            use_vaxijen: 是否使用完整VaxiJen算法（推荐True）
        """
        self.protegenicity_threshold = protegenicity_threshold
        self.antigenicity_threshold = antigenicity_threshold
        self.organism_type = organism_type
        self.use_vaxijen = use_vaxijen
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # 初始化VaxiJen计算器
        if self.use_vaxijen:
            self.vaxijen_calculator = VaxiJenCalculator(
                organism_type=organism_type,
                custom_threshold=antigenicity_threshold
            )
            self.logger.info(f"使用VaxiJen 2.0算法 (生物类型: {organism_type})")
    
    def evaluate_from_sequences(self,
                                sequences: Dict[str, ProteinSequence]) -> pd.DataFrame:
        """
        基于VaxiJen 2.0算法评估蛋白质免疫原性（不依赖Vaxign-ML）

        对每个蛋白质使用VaxiJen 2.0计算抗原性分数，
        将其映射为0-100的保护性分数用于候选排序。

        Args:
            sequences: 序列字典 {protein_id: ProteinSequence}

        Returns:
            免疫原性评估结果DataFrame
        """
        self.logger.info(f"开始免疫原性评估 (VaxiJen 2.0)，{len(sequences)} 个蛋白质")

        results = []

        for protein_id, sequence in sequences.items():
            # 使用VaxiJen 2.0计算抗原性分数
            antigenicity = self._calculate_antigenicity(sequence)

            # 将抗原性分数(0-1)映射为保护性分数(0-100)
            protegenicity = antigenicity * 100

            # 判断是否为可能的保护性抗原
            is_protective = protegenicity >= self.protegenicity_threshold

            result = ImmunogenicityResult(
                protein_id=protein_id,
                protegenicity_score=protegenicity,
                is_protective_antigen=is_protective,
                subcellular_location='Unknown',
                has_signal_peptide=False,
                transmembrane_regions=0,
                antigenicity_score=antigenicity,
                immunogenicity_rank=0,
                recommendation='',
                features={'antigenicity_score': antigenicity}
            )

            results.append(result)

        # 计算排名和推荐
        results = self._rank_and_recommend(results)

        # 转换为DataFrame
        df = self._results_to_dataframe(results)

        self.logger.info(f"免疫原性评估完成，{len(df)} 个蛋白质")

        return df

    def evaluate(self,
                 vaxignml_result,
                 sequences: Dict[str, ProteinSequence]) -> pd.DataFrame:
        """
        评估蛋白质的免疫原性（兼容旧接口，需Vaxign-ML结果）

        Args:
            vaxignml_result: Vaxign-ML预测结果
            sequences: 序列字典

        Returns:
            免疫原性评估结果DataFrame
        """
        self.logger.info("开始免疫原性评估")

        results = []

        # 从预测结果中提取蛋白质级别分数
        protein_scores = vaxignml_result.protein_scores

        for protein_id, protegenicity in protein_scores.items():
            # 从预测元数据中提取额外信息
            features = self._extract_features(
                protein_id,
                vaxignml_result.predictions,
                sequences.get(protein_id)
            )

            # 计算抗原性分数
            antigenicity = self._calculate_antigenicity(
                sequences.get(protein_id)
            )

            # 确定亚细胞定位
            location = features.get('subcellular_location', 'Unknown')

            # 判断是否为保护性抗原
            is_protective = protegenicity >= self.protegenicity_threshold

            result = ImmunogenicityResult(
                protein_id=protein_id,
                protegenicity_score=protegenicity,
                is_protective_antigen=is_protective,
                subcellular_location=location,
                has_signal_peptide=features.get('has_signal_peptide', False),
                transmembrane_regions=features.get('transmembrane_regions', 0),
                antigenicity_score=antigenicity,
                immunogenicity_rank=0,  # 后续计算
                recommendation='',  # 后续确定
                features=features
            )

            results.append(result)

        # 计算排名和推荐
        results = self._rank_and_recommend(results)

        # 转换为DataFrame
        df = self._results_to_dataframe(results)

        self.logger.info(f"免疫原性评估完成，{len(df)} 个蛋白质")

        return df
    
    def _extract_features(self, 
                           protein_id: str,
                           predictions: list,
                           sequence: Optional[ProteinSequence]) -> dict:
        """
        从预测结果中提取特征
        """
        features = {}
        
        # 查找该蛋白质的预测
        protein_preds = [p for p in predictions if p.protein_id == protein_id]
        
        if protein_preds:
            # 从第一个预测中提取额外信息
            first_pred = protein_preds[0]
            additional = first_pred.additional_info
            
            features['subcellular_location'] = additional.get(
                'subcellular_location', 
                'Unknown'
            )
            features['has_signal_peptide'] = additional.get(
                'has_signal_peptide', 
                False
            )
            features['transmembrane_regions'] = additional.get(
                'transmembrane_regions', 
                0
            )
            features['protegenicity'] = additional.get('protegenicity', 0)
        
        return features
    
    def _calculate_antigenicity(self, 
                                 sequence: Optional[ProteinSequence]) -> float:
        """
        计算抗原性分数
        
        使用完整VaxiJen 2.0算法或简化方法
        """
        if not sequence:
            return 0.0
        
        seq = sequence.sequence.upper()
        n = len(seq)
        
        if n == 0:
            return 0.0
        
        # 使用完整VaxiJen算法
        if self.use_vaxijen and hasattr(self, 'vaxijen_calculator'):
            vaxijen_result = self.vaxijen_calculator.calculate(sequence)
            return vaxijen_result.antigenicity_score
        
        # 备用: 简化的氨基酸组成计算
        return self._calculate_simple_antigenicity(sequence)
    
    def _calculate_simple_antigenicity(self, 
                                       sequence: ProteinSequence) -> float:
        """
        简化的抗原性计算（备用方法）
        
        使用Welling抗原性规模
        """
        seq = sequence.sequence.upper()
        
        # Welling抗原性系数
        antigenicity_coefs = {
            'A': 1.064, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091,
            'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 1.250,
            'M': 0.826, 'N': 0.776, 'P': 1.064, 'Q': 1.015, 'R': 0.873,
            'S': 0.883, 'T': 0.909, 'V': 1.383, 'W': 0.893, 'Y': 1.161
        }
        
        total = 0.0
        count = 0
        
        for aa in seq:
            if aa in antigenicity_coefs:
                total += antigenicity_coefs[aa]
                count += 1
        
        if count == 0:
            return 0.0
        
        # 归一化到0-1范围
        avg = total / count
        normalized = (avg - 0.7) / 0.7
        return max(0.0, min(1.0, normalized))
    
    def get_vaxijen_details(self, 
                            sequence: ProteinSequence) -> Optional[Dict]:
        """
        获取VaxiJen详细分析结果
        
        Args:
            sequence: 蛋白质序列
        
        Returns:
            详细分析字典，如果未启用VaxiJen则返回None
        """
        if self.use_vaxijen and hasattr(self, 'vaxijen_calculator'):
            return self.vaxijen_calculator.get_detailed_analysis(sequence)
        return None
    
    def _rank_and_recommend(self,
                            results: List[ImmunogenicityResult]) -> List[ImmunogenicityResult]:
        """
        计算排名和推荐等级

        基于VaxiJen抗原性分数和保护性分数进行排序和推荐。
        """
        if not results:
            return results

        # 计算综合评分
        for result in results:
            # 综合评分 = 保护性分数 * 抗原性分数
            # 保护性分数已由VaxiJen抗原性映射(0-100)，抗原性分数为0-1
            result.features['composite_score'] = (
                result.protegenicity_score *
                result.antigenicity_score
            )

        # 按综合评分排序
        results.sort(key=lambda x: x.features.get('composite_score', 0), reverse=True)

        # 分配排名和推荐
        for i, result in enumerate(results):
            result.immunogenicity_rank = i + 1

            # 确定推荐等级（基于VaxiJen抗原性分数）
            if result.antigenicity_score >= 0.7 and result.is_protective_antigen:
                result.recommendation = 'HIGH'
            elif result.antigenicity_score >= self.antigenicity_threshold:
                result.recommendation = 'MEDIUM'
            else:
                result.recommendation = 'LOW'

        return results
    
    def _get_location_weight(self, location: str) -> float:
        """获取亚细胞定位权重"""
        weights = {
            'OuterMembrane': 1.5,
            'Extracellular': 1.4,
            'Periplasmic': 1.3,
            'CellWall': 1.3,
            'Secreted': 1.4,
            'Fimbrial': 1.5,
            'Surface': 1.5,
            'InnerMembrane': 1.0,
            'Cytoplasmic': 0.8,
            'Unknown': 1.0
        }
        
        return weights.get(location, 1.0)
    
    def _results_to_dataframe(self, 
                               results: List[ImmunogenicityResult]) -> pd.DataFrame:
        """将结果转换为DataFrame"""
        data = []
        
        for result in results:
            row = {
                'protein_id': result.protein_id,
                'protegenicity_score': result.protegenicity_score,
                'is_protective_antigen': result.is_protective_antigen,
                'subcellular_location': result.subcellular_location,
                'has_signal_peptide': result.has_signal_peptide,
                'transmembrane_regions': result.transmembrane_regions,
                'antigenicity_score': result.antigenicity_score,
                'immunogenicity_rank': result.immunogenicity_rank,
                'recommendation': result.recommendation,
                'composite_score': result.features.get('composite_score', 0)
            }
            data.append(row)
        
        return pd.DataFrame(data)
    
    def get_top_candidates(self, 
                           immunogenicity_df: pd.DataFrame,
                           n: int = 50,
                           recommendation: str = None) -> pd.DataFrame:
        """
        获取Top N候选抗原
        
        Args:
            immunogenicity_df: 免疫原性结果DataFrame
            n: 返回数量
            recommendation: 筛选推荐等级（可选）
        
        Returns:
            Top N候选DataFrame
        """
        df = immunogenicity_df.copy()
        
        if recommendation:
            df = df[df['recommendation'] == recommendation]
        
        return df.head(n)
