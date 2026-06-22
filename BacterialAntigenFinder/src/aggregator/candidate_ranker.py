"""
候选抗原排序模块

整合共识评分和免疫原性评估结果，生成最终的候选抗原排名。
"""

import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
import pandas as pd
import numpy as np

from ..preprocessor.fasta_parser import ProteinSequence


@dataclass
class AntigenCandidate:
    """抗原候选数据类"""
    protein_id: str
    epitope_regions: List[Tuple[int, int, float, str]]  # (start, end, score, sequence)
    total_epitope_residues: int
    avg_consensus_score: float
    max_consensus_score: float
    protegenicity_score: float
    immunogenicity_rank: int
    subcellular_location: str
    recommendation: str
    composite_score: float
    predictor_agreement: Dict[str, int] = field(default_factory=dict)


class CandidateRanker:
    """候选抗原排序器"""
    
    def __init__(self, 
                 top_n: int = 50,
                 min_epitope_length: int = 5,
                 min_consensus_score: float = 0.3):
        """
        初始化排序器
        
        Args:
            top_n: 输出前N个候选
            min_epitope_length: 最小表位区域长度
            min_consensus_score: 最小共识分数
        """
        self.top_n = top_n
        self.min_epitope_length = min_epitope_length
        self.min_consensus_score = min_consensus_score
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def rank(self, 
             consensus_df: pd.DataFrame,
             immunogenicity_df: Optional[pd.DataFrame],
             sequences: Dict[str, ProteinSequence]) -> pd.DataFrame:
        """
        对候选抗原进行排序
        
        Args:
            consensus_df: 共识评分结果
            immunogenicity_df: 免疫原性评估结果
            sequences: 序列字典
        
        Returns:
            排序后的候选抗原DataFrame
        """
        self.logger.info("开始候选抗原排序")
        
        # 空数据检查
        if consensus_df is None or consensus_df.empty:
            self.logger.warning("共识评分为空，跳过候选排序")
            return pd.DataFrame()
        
        candidates = []
        
        # 获取所有蛋白质ID
        protein_ids = consensus_df['protein_id'].unique()
        
        for protein_id in protein_ids:
            # 提取该蛋白质的共识结果
            protein_consensus = consensus_df[
                consensus_df['protein_id'] == protein_id
            ]
            
            # 提取表位区域
            epitope_regions = self._extract_epitope_regions(protein_consensus)
            
            if not epitope_regions:
                continue
            
            # 计算统计信息
            epitope_residues = protein_consensus[
                protein_consensus['is_consensus_epitope'] == True
            ]
            
            total_epitope_residues = len(epitope_residues)
            avg_score = epitope_residues['consensus_score'].mean() if not epitope_residues.empty else 0
            max_score = epitope_residues['consensus_score'].max() if not epitope_residues.empty else 0
            
            # 获取免疫原性信息
            if immunogenicity_df is not None and not immunogenicity_df.empty:
                immuno_row = immunogenicity_df[
                    immunogenicity_df['protein_id'] == protein_id
                ]
                if not immuno_row.empty:
                    protegenicity = immuno_row['protegenicity_score'].values[0]
                    immuno_rank = immuno_row['immunogenicity_rank'].values[0]
                    location = immuno_row['subcellular_location'].values[0]
                    recommendation = immuno_row['recommendation'].values[0]
                else:
                    protegenicity = 0
                    immuno_rank = 999
                    location = 'Unknown'
                    recommendation = 'LOW'
            else:
                protegenicity = 0
                immuno_rank = 999
                location = 'Unknown'
                recommendation = 'LOW'
            
            # 计算预测器一致性
            predictor_cols = [col for col in protein_consensus.columns 
                             if col.endswith('_epitope')]
            predictor_agreement = {}
            for col in predictor_cols:
                predictor_name = col.replace('_epitope', '')
                votes = protein_consensus[col].sum()
                predictor_agreement[predictor_name] = int(votes)
            
            # 计算综合评分
            composite_score = self._calculate_composite_score(
                avg_score, max_score, total_epitope_residues,
                protegenicity, len(epitope_regions)
            )
            
            candidate = AntigenCandidate(
                protein_id=protein_id,
                epitope_regions=epitope_regions,
                total_epitope_residues=total_epitope_residues,
                avg_consensus_score=avg_score,
                max_consensus_score=max_score,
                protegenicity_score=protegenicity,
                immunogenicity_rank=immuno_rank,
                subcellular_location=location,
                recommendation=recommendation,
                composite_score=composite_score,
                predictor_agreement=predictor_agreement
            )
            
            candidates.append(candidate)
        
        # 按综合评分排序
        candidates.sort(key=lambda x: x.composite_score, reverse=True)
        
        # 转换为DataFrame
        df = self._candidates_to_dataframe(candidates, sequences)
        
        # 返回Top N
        result_df = df.head(self.top_n)
        
        self.logger.info(f"候选排序完成，返回 {len(result_df)} 个候选")
        
        return result_df
    
    def _extract_epitope_regions(self, 
                                  protein_df: pd.DataFrame) -> List[Tuple[int, int, float, str]]:
        """提取连续的表位区域"""
        epitope_df = protein_df[
            protein_df['is_consensus_epitope'] == True
        ].sort_values('residue_id')
        
        if epitope_df.empty:
            return []
        
        regions = []
        current_start = None
        current_end = None
        current_scores = []
        current_seq = ""
        
        for _, row in epitope_df.iterrows():
            res_id = row['residue_id']
            
            if current_start is None:
                current_start = res_id
                current_end = res_id
                current_scores = [row['consensus_score']]
                current_seq = row['residue_name']
            elif res_id == current_end + 1:
                current_end = res_id
                current_scores.append(row['consensus_score'])
                current_seq += row['residue_name']
            else:
                # 保存当前区域
                if current_end - current_start + 1 >= self.min_epitope_length:
                    avg_score = sum(current_scores) / len(current_scores)
                    if avg_score >= self.min_consensus_score:
                        regions.append((current_start, current_end, avg_score, current_seq))
                
                # 开始新区域
                current_start = res_id
                current_end = res_id
                current_scores = [row['consensus_score']]
                current_seq = row['residue_name']
        
        # 处理最后一个区域
        if current_start is not None and current_end - current_start + 1 >= self.min_epitope_length:
            avg_score = sum(current_scores) / len(current_scores)
            if avg_score >= self.min_consensus_score:
                regions.append((current_start, current_end, avg_score, current_seq))
        
        return regions
    
    def _calculate_composite_score(self,
                                    avg_score: float,
                                    max_score: float,
                                    total_residues: int,
                                    protegenicity: float,
                                    num_regions: int) -> float:
        """计算综合评分"""
        # 归一化各分量
        score_component = (avg_score + max_score) / 2
        
        # 残基数量分量（对数缩放）
        residue_component = np.log1p(total_residues) / 5  # 归一化到大约0-1
        
        # 保护性分数分量
        protegen_component = protegenicity / 100 if protegenicity else 0
        
        # 区域数量分量
        region_component = min(num_regions / 5, 1.0)  # 最多考虑5个区域
        
        # 加权综合
        composite = (
            0.35 * score_component +
            0.25 * protegen_component +
            0.20 * residue_component +
            0.20 * region_component
        )
        
        return composite
    
    def _candidates_to_dataframe(self, 
                                  candidates: List[AntigenCandidate],
                                  sequences: Dict[str, ProteinSequence]) -> pd.DataFrame:
        """将候选列表转换为DataFrame"""
        data = []
        
        for i, candidate in enumerate(candidates):
            # 获取最佳表位区域
            if candidate.epitope_regions:
                best_region = max(candidate.epitope_regions, key=lambda x: x[2])
                best_start, best_end, best_score, best_seq = best_region
                residue_range = f"{best_start}-{best_end}"
                epitope_sequence = best_seq
            else:
                residue_range = ""
                epitope_sequence = ""
            
            # 获取序列信息
            seq_obj = sequences.get(candidate.protein_id)
            sequence_length = seq_obj.length if seq_obj else 0
            
            row = {
                'rank': i + 1,
                'protein_id': candidate.protein_id,
                'sequence_length': sequence_length,
                'residue_range': residue_range,
                'epitope_sequence': epitope_sequence,
                'num_epitope_regions': len(candidate.epitope_regions),
                'total_epitope_residues': candidate.total_epitope_residues,
                'avg_consensus_score': round(candidate.avg_consensus_score, 4),
                'max_consensus_score': round(candidate.max_consensus_score, 4),
                'protegenicity_score': round(candidate.protegenicity_score, 2),
                'immunogenicity_rank': candidate.immunogenicity_rank,
                'subcellular_location': candidate.subcellular_location,
                'recommendation': candidate.recommendation,
                'composite_score': round(candidate.composite_score, 4)
            }
            
            # 添加各预测器投票数
            for predictor, votes in candidate.predictor_agreement.items():
                row[f'{predictor}_votes'] = votes
            
            # 添加所有表位区域
            all_regions = '; '.join([
                f"{s}-{e}({score:.2f})" 
                for s, e, score, seq in candidate.epitope_regions
            ])
            row['all_epitope_regions'] = all_regions
            
            data.append(row)
        
        return pd.DataFrame(data)
    
    def rank_epitopes(self,
                      consensus_df: pd.DataFrame,
                      immunogenicity_df: Optional[pd.DataFrame],
                      sequences: Dict[str, ProteinSequence]) -> pd.DataFrame:
        """
        对表位级别进行排序

        与rank()不同，此方法返回每一行是一个独立的表位区域，
        而非蛋白质级别的汇总。

        Args:
            consensus_df: 共识评分结果
            immunogenicity_df: 免疫原性评估结果
            sequences: 序列字典

        Returns:
            排序后的表位级别DataFrame
        """
        self.logger.info("开始表位级别排序")

        # 空数据检查
        if consensus_df is None or consensus_df.empty:
            self.logger.warning("共识评分为空，跳过表位排序")
            return pd.DataFrame()

        epitope_rows = []

        # 获取所有蛋白质ID
        protein_ids = consensus_df['protein_id'].unique()

        # 识别预测器列
        predictor_score_cols = [col for col in consensus_df.columns
                                if col.endswith('_score') and col != 'consensus_score']
        predictor_epitope_cols = [col for col in consensus_df.columns
                                  if col.endswith('_epitope') and col != 'is_consensus_epitope']
        predictor_names = [col.replace('_score', '') for col in predictor_score_cols]

        for protein_id in protein_ids:
            # 提取该蛋白质的共识结果
            protein_consensus = consensus_df[
                consensus_df['protein_id'] == protein_id
            ]

            # 提取表位区域
            epitope_regions = self._extract_epitope_regions(protein_consensus)

            if not epitope_regions:
                continue

            # 获取免疫原性信息
            if immunogenicity_df is not None and not immunogenicity_df.empty:
                immuno_row = immunogenicity_df[
                    immunogenicity_df['protein_id'] == protein_id
                ]
                if not immuno_row.empty:
                    protegenicity = immuno_row['protegenicity_score'].values[0]
                    immuno_rank = immuno_row['immunogenicity_rank'].values[0]
                    location = immuno_row['subcellular_location'].values[0]
                    recommendation = immuno_row['recommendation'].values[0]
                    antigenicity = immuno_row['antigenicity_score'].values[0] \
                        if 'antigenicity_score' in immuno_row.columns else 0.0
                else:
                    protegenicity = 0
                    immuno_rank = 999
                    location = 'Unknown'
                    recommendation = 'LOW'
                    antigenicity = 0.0
            else:
                protegenicity = 0
                immuno_rank = 999
                location = 'Unknown'
                recommendation = 'LOW'
                antigenicity = 0.0

            for start, end, avg_score, epitope_seq in epitope_regions:
                # 提取该表位区域内的残基数据
                region_residues = protein_consensus[
                    (protein_consensus['residue_id'] >= start) &
                    (protein_consensus['residue_id'] <= end)
                ]

                epitope_length = end - start + 1
                max_consensus = region_residues['consensus_score'].max() \
                    if not region_residues.empty else 0
                total_votes = int(region_residues['vote_count'].sum()) \
                    if not region_residues.empty else 0

                # 计算预测器一致性：在该表位区域内，有多少预测器在至少一个残基上判定为表位
                num_predictors_agree = 0
                predictor_scores = {}
                for pname in predictor_names:
                    epitope_col = f'{pname}_epitope'
                    score_col = f'{pname}_score'
                    if epitope_col in region_residues.columns:
                        if region_residues[epitope_col].any():
                            num_predictors_agree += 1
                    # 收集该预测器在表位区域内的平均分数
                    if score_col in region_residues.columns:
                        predictor_scores[f'{pname}_score'] = \
                            region_residues[score_col].mean()
                    else:
                        predictor_scores[f'{pname}_score'] = None

                # 计算表位级别的综合评分
                composite_score = self._calculate_epitope_composite_score(
                    avg_consensus_score=avg_score,
                    protegenicity=protegenicity,
                    antigenicity=antigenicity,
                    epitope_length=epitope_length,
                    num_predictors_agree=num_predictors_agree,
                    total_predictors=len(predictor_names) if predictor_names else 1
                )

                # 确定推荐等级
                epitope_recommendation = self._determine_epitope_recommendation(
                    composite_score, avg_score, protegenicity, antigenicity
                )

                row = {
                    'protein_id': protein_id,
                    'epitope_start': start,
                    'epitope_end': end,
                    'epitope_sequence': epitope_seq,
                    'epitope_length': epitope_length,
                    'avg_consensus_score': round(avg_score, 4),
                    'max_consensus_score': round(max_consensus, 4),
                    'vote_count': total_votes,
                    'num_predictors_agree': num_predictors_agree,
                    'protegenicity_score': round(protegenicity, 2),
                    'antigenicity_score': round(antigenicity, 4),
                    'subcellular_location': location,
                    'immunogenicity_rank': immuno_rank,
                    'recommendation': epitope_recommendation,
                    'composite_score': round(composite_score, 4),
                }

                # 添加各预测器分数
                row.update(predictor_scores)

                epitope_rows.append(row)

        if not epitope_rows:
            self.logger.info("未找到符合条件的表位")
            return pd.DataFrame()

        # 构建DataFrame并排序
        result_df = pd.DataFrame(epitope_rows)
        result_df.sort_values('composite_score', ascending=False, inplace=True)
        result_df.reset_index(drop=True, inplace=True)
        result_df['rank'] = range(1, len(result_df) + 1)

        # 返回Top N
        result_df = result_df.head(self.top_n)

        self.logger.info(f"表位级别排序完成，返回 {len(result_df)} 个表位")

        return result_df

    def _calculate_epitope_composite_score(self,
                                           avg_consensus_score: float,
                                           protegenicity: float,
                                           antigenicity: float,
                                           epitope_length: int,
                                           num_predictors_agree: int,
                                           total_predictors: int) -> float:
        """
        计算表位级别的综合评分

        权重分配:
        - 0.30 * normalized_avg_consensus_score
        - 0.25 * normalized_protegenicity
        - 0.20 * normalized_antigenicity
        - 0.15 * epitope_length_factor (对数缩放，上限1.0)
        - 0.10 * predictor_agreement_factor

        Args:
            avg_consensus_score: 表位区域平均共识分数
            protegenicity: 保护性分数
            antigenicity: 抗原性分数
            epitope_length: 表位长度
            num_predictors_agree: 一致同意的预测器数量
            total_predictors: 总预测器数量

        Returns:
            综合评分
        """
        # 归一化共识分数（通常在0-1之间）
        normalized_consensus = min(max(avg_consensus_score, 0.0), 1.0)

        # 归一化保护性分数（通常0-100）
        normalized_protegenicity = min(protegenicity / 100.0, 1.0) if protegenicity else 0.0

        # 归一化抗原性分数（通常0-1）
        normalized_antigenicity = min(max(antigenicity, 0.0), 1.0)

        # 表位长度因子：对数缩放，上限1.0
        epitope_length_factor = min(np.log1p(epitope_length) / np.log1p(30), 1.0)

        # 预测器一致性因子
        predictor_agreement_factor = num_predictors_agree / total_predictors \
            if total_predictors > 0 else 0.0

        composite = (
            0.30 * normalized_consensus +
            0.25 * normalized_protegenicity +
            0.20 * normalized_antigenicity +
            0.15 * epitope_length_factor +
            0.10 * predictor_agreement_factor
        )

        return composite

    def _determine_epitope_recommendation(self,
                                          composite_score: float,
                                          avg_consensus_score: float,
                                          protegenicity: float,
                                          antigenicity: float) -> str:
        """
        根据综合评分和各分量确定表位推荐等级

        Args:
            composite_score: 综合评分
            avg_consensus_score: 平均共识分数
            protegenicity: 保护性分数
            antigenicity: 抗原性分数

        Returns:
            推荐等级: HIGH / MEDIUM / LOW
        """
        if composite_score >= 0.6 and avg_consensus_score >= 0.5 and protegenicity >= 50:
            return 'HIGH'
        elif composite_score >= 0.4 and (avg_consensus_score >= 0.3 or protegenicity >= 30):
            return 'MEDIUM'
        else:
            return 'LOW'

    def get_epitope_details(self,
                           candidates_df: pd.DataFrame,
                           protein_id: str) -> Dict:
        """获取特定蛋白质的表位详情"""
        row = candidates_df[candidates_df['protein_id'] == protein_id]
        
        if row.empty:
            return {}
        
        row = row.iloc[0]
        
        return {
            'protein_id': protein_id,
            'rank': row['rank'],
            'epitope_regions': row['all_epitope_regions'],
            'total_residues': row['total_epitope_residues'],
            'avg_score': row['avg_consensus_score'],
            'protegenicity': row['protegenicity_score'],
            'location': row['subcellular_location'],
            'recommendation': row['recommendation']
        }
