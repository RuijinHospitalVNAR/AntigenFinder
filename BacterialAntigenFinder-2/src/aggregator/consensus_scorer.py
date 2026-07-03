"""
共识评分引擎

整合多个预测器的结果，计算加权共识分数。
"""

import logging
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from dataclasses import dataclass, field
import pandas as pd
import numpy as np

from ..predictors.base_predictor import PredictionResult, EpitopePrediction


@dataclass
class ConsensusResult:
    """共识评分结果"""
    protein_id: str
    residue_id: int
    residue_name: str
    consensus_score: float
    vote_count: int
    is_consensus_epitope: bool
    individual_scores: Dict[str, float] = field(default_factory=dict)
    individual_predictions: Dict[str, bool] = field(default_factory=dict)


class ConsensusScorer:
    """共识评分器
    
    支持两种分数归一化模式:
    - 'none': 使用原始分数 (仅当所有模型输出相同尺度时适用)
    - 'rank': 每个模型内部转换为百分位排名 [0,1]，使不同模型的分数可比
    """
    
    DEFAULT_WEIGHTS = {
        'bepipred': 0.30,
        'discotope': 0.30,
        'graphbepi': 0.20,
        'epigraph': 0.20
    }
    
    def __init__(self, 
                 weights: Optional[Dict[str, float]] = None,
                 threshold: float = 0.5,
                 min_votes: int = 2,
                 method: str = 'weighted_avg',
                 score_normalize: str = 'rank'):
        """
        初始化共识评分器
        
        Args:
            weights: 各预测器权重
            threshold: 共识评分阈值 (归一化后 [0,1] 区间)
            min_votes: 最小投票数
            method: 评分方法 ('weighted_avg', 'majority_vote', 'max_score', 'mean')
            score_normalize: 分数归一化方法 ('none' | 'rank')
        """
        self.weights = weights or self.DEFAULT_WEIGHTS
        self.threshold = threshold
        self.min_votes = min_votes
        self.method = method
        self.score_normalize = score_normalize
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def compute_consensus(self, 
                          predictions: Dict[str, PredictionResult]) -> pd.DataFrame:
        """
        计算共识评分
        
        Args:
            predictions: 预测器名称到预测结果的字典
        
        Returns:
            共识结果DataFrame
        """
        if not predictions:
            self.logger.warning("没有预测结果")
            return pd.DataFrame()
        
        # 自适应调整min_votes：不超过实际预测器数量
        num_predictors = len(predictions)
        effective_min_votes = min(self.min_votes, num_predictors)
        if effective_min_votes < self.min_votes:
            self.logger.info(
                f"预测器数量({num_predictors})少于min_votes({self.min_votes})，"
                f"自动调整为{effective_min_votes}"
            )
        
        self.logger.info(f"开始计算共识评分，{num_predictors} 个预测器，min_votes={effective_min_votes}")
        
        # 收集所有残基的预测
        residue_predictions = self._collect_residue_predictions(predictions)
        
        # 分数归一化: 将每个模型的原始分数转为百分位排名
        if self.score_normalize == 'rank':
            residue_predictions = self._normalize_scores_by_rank(residue_predictions)
            self.logger.info("已对各模型分数进行百分位排名归一化")
        
        # 计算共识分数
        consensus_results = []
        
        for (protein_id, residue_id), pred_dict in residue_predictions.items():
            result = self._compute_single_consensus(
                protein_id, residue_id, pred_dict, effective_min_votes
            )
            if result:
                consensus_results.append(result)
        
        # 转换为DataFrame
        df = self._results_to_dataframe(consensus_results)
        
        self.logger.info(f"共识评分完成，{len(df)} 个残基")
        
        return df
    
    def _collect_residue_predictions(self, 
                                      predictions: Dict[str, PredictionResult]) -> dict:
        """
        收集每个残基的所有预测
        
        Returns:
            {(protein_id, residue_id): {predictor: (score, is_epitope, residue_name)}}
        """
        residue_predictions = defaultdict(dict)
        
        for predictor_name, pred_result in predictions.items():
            for pred in pred_result.predictions:
                key = (pred.protein_id, pred.residue_id)
                residue_predictions[key][predictor_name] = (
                    pred.score,
                    pred.is_epitope,
                    pred.residue_name
                )
        
        return residue_predictions
    
    def _normalize_scores_by_rank(self, 
                                   residue_predictions: dict) -> dict:
        """
        将每个预测器的所有残基分数转换为百分位排名 [0,1]。
        
        变换后各模型分数在相同尺度上可比，解决了 BepiPred(0-1概率)、
        GraphBepi(GNN logit 0-0.4)、EpiGraph(GAT logit 0-0.25) 等
        分数尺度不一致的问题。
        """
        # 1. 收集每个预测器的所有原始分数
        predictor_raw_scores = defaultdict(list)
        for (protein_id, residue_id), pred_dict in residue_predictions.items():
            for predictor, (score, is_epitope, res_name) in pred_dict.items():
                predictor_raw_scores[predictor].append(score)
        
        # 2. 计算百分位排名映射
        predictor_rank_map = {}
        for predictor, scores in predictor_raw_scores.items():
            arr = np.array(scores)
            n = len(arr)
            if n <= 1:
                # 只有一个残基时无法计算排名，赋值为 0.5
                for key in residue_predictions:
                    if predictor in residue_predictions[key]:
                        predictor_rank_map[(key[0], key[1], predictor)] = 0.5
                continue
            
            # 纯 NumPy 实现 rankdata(arr, method='average')
            # 排序并计算每个元素的排名（平均值处理并列）
            order = np.argsort(arr)
            ranks = np.empty(n, dtype=float)
            i = 0
            while i < n:
                j = i
                # 找到并列值的范围
                while j < n and arr[order[j]] == arr[order[i]]:
                    j += 1
                # 并列值使用平均排名
                avg_rank = (i + j - 1) / 2.0  # 0-based average rank of the group
                for k in range(i, j):
                    ranks[order[k]] = avg_rank
                i = j
            # 转为 [0, 1] 百分位
            percentiles = ranks / (n - 1) if n > 1 else np.zeros(n)
            # 建立 (protein_id, residue_id) -> percentile 的快速查找
            idx = 0
            for key in residue_predictions:
                if predictor in residue_predictions[key]:
                    predictor_rank_map[(key[0], key[1], predictor)] = float(percentiles[idx])
                    idx += 1
            self.logger.debug(
                f"{predictor}: 分数范围 [{min(scores):.4f}, {max(scores):.4f}] -> "
                f"排名范围 [0, 1], n={n}"
            )
        
        # 3. 替换分数为百分位排名
        normalized = {}
        for (protein_id, residue_id), pred_dict in residue_predictions.items():
            new_dict = {}
            for predictor, (score, is_epitope, res_name) in pred_dict.items():
                norm_score = predictor_rank_map.get(
                    (protein_id, residue_id, predictor), 0.0
                )
                new_dict[predictor] = (norm_score, is_epitope, res_name)
            normalized[(protein_id, residue_id)] = new_dict
        
        return normalized
    
    def _compute_single_consensus(self, 
                                   protein_id: str,
                                   residue_id: int,
                                   pred_dict: dict,
                                   effective_min_votes: int = None) -> Optional[ConsensusResult]:
        """
        计算单个残基的共识分数
        """
        if not pred_dict:
            return None
        
        # 使用传入的min_votes或默认值
        min_votes = effective_min_votes if effective_min_votes is not None else self.min_votes
        
        # 收集分数和投票
        scores = {}
        votes = {}
        residue_name = 'X'
        
        for predictor, (score, is_epitope, res_name) in pred_dict.items():
            scores[predictor] = score
            votes[predictor] = is_epitope
            if res_name and res_name != 'X':
                residue_name = res_name
        
        # 计算共识分数
        if self.method == 'weighted_avg':
            consensus_score = self._weighted_average(scores)
        elif self.method == 'majority_vote':
            consensus_score = self._majority_vote_score(votes)
        elif self.method == 'max_score':
            consensus_score = max(scores.values()) if scores else 0.0
        else:
            consensus_score = self._weighted_average(scores)
        
        # 计算投票数
        vote_count = sum(1 for v in votes.values() if v)
        
        # 判断是否为共识表位
        is_consensus_epitope = (
            consensus_score >= self.threshold and 
            vote_count >= min_votes
        )
        
        return ConsensusResult(
            protein_id=protein_id,
            residue_id=residue_id,
            residue_name=residue_name,
            consensus_score=consensus_score,
            vote_count=vote_count,
            is_consensus_epitope=is_consensus_epitope,
            individual_scores=scores,
            individual_predictions=votes
        )
    
    def _weighted_average(self, scores: Dict[str, float]) -> float:
        """计算加权平均分数"""
        if not scores:
            return 0.0
        
        total_weight = 0.0
        weighted_sum = 0.0
        
        for predictor, score in scores.items():
            weight = self.weights.get(predictor, 0.1)
            weighted_sum += score * weight
            total_weight += weight
        
        if total_weight == 0:
            return sum(scores.values()) / len(scores)
        
        return weighted_sum / total_weight
    
    def _majority_vote_score(self, votes: Dict[str, bool]) -> float:
        """计算多数投票分数"""
        if not votes:
            return 0.0
        
        vote_count = sum(1 for v in votes.values() if v)
        return vote_count / len(votes)
    
    def _results_to_dataframe(self, results: List[ConsensusResult]) -> pd.DataFrame:
        """将结果转换为DataFrame"""
        data = []
        
        for result in results:
            row = {
                'protein_id': result.protein_id,
                'residue_id': result.residue_id,
                'residue_name': result.residue_name,
                'consensus_score': result.consensus_score,
                'vote_count': result.vote_count,
                'is_consensus_epitope': result.is_consensus_epitope
            }
            
            # 添加各预测器分数
            for predictor in self.weights.keys():
                row[f'{predictor}_score'] = result.individual_scores.get(predictor, None)
                row[f'{predictor}_epitope'] = result.individual_predictions.get(predictor, None)
            
            data.append(row)
        
        df = pd.DataFrame(data)
        
        # 按蛋白质ID和残基ID排序
        if not df.empty:
            df = df.sort_values(['protein_id', 'residue_id'])
        
        return df
    
    def get_epitope_regions(self, 
                            consensus_df: pd.DataFrame,
                            protein_id: str,
                            min_length: int = 5) -> List[Tuple[int, int, float, str]]:
        """
        获取连续的表位区域
        
        Args:
            consensus_df: 共识结果DataFrame
            protein_id: 蛋白质ID
            min_length: 最小区域长度
        
        Returns:
            [(start, end, avg_score, sequence), ...]
        """
        # 筛选该蛋白质的表位
        protein_df = consensus_df[
            (consensus_df['protein_id'] == protein_id) & 
            (consensus_df['is_consensus_epitope'] == True)
        ].sort_values('residue_id')
        
        if protein_df.empty:
            return []
        
        regions = []
        current_start = None
        current_end = None
        current_scores = []
        current_seq = ""
        
        for _, row in protein_df.iterrows():
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
                if current_end - current_start + 1 >= min_length:
                    avg_score = sum(current_scores) / len(current_scores)
                    regions.append((current_start, current_end, avg_score, current_seq))
                
                # 开始新区域
                current_start = res_id
                current_end = res_id
                current_scores = [row['consensus_score']]
                current_seq = row['residue_name']
        
        # 处理最后一个区域
        if current_start is not None and current_end - current_start + 1 >= min_length:
            avg_score = sum(current_scores) / len(current_scores)
            regions.append((current_start, current_end, avg_score, current_seq))
        
        return regions
    
    def get_summary_stats(self, consensus_df: pd.DataFrame) -> dict:
        """获取共识结果的统计摘要"""
        if consensus_df.empty:
            return {}
        
        return {
            'total_residues': len(consensus_df),
            'epitope_residues': consensus_df['is_consensus_epitope'].sum(),
            'epitope_ratio': consensus_df['is_consensus_epitope'].mean(),
            'avg_consensus_score': consensus_df['consensus_score'].mean(),
            'max_consensus_score': consensus_df['consensus_score'].max(),
            'avg_vote_count': consensus_df['vote_count'].mean(),
            'proteins': consensus_df['protein_id'].nunique()
        }
