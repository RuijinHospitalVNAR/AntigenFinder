"""
BacterialAntigenFinder - 细菌抗原AI智能筛选平台

整合多种B细胞表位预测模型，通过共识评分和免疫原性评估，
输出针对耐药细菌的抗原候选清单。

支持的模型:
- BepiPred-3.0: 基于序列的B细胞表位预测
- DiscoTope-3.0: 基于结构的构象表位预测
- GraphBepi: 图神经网络表位预测
- EpiGraph: 图注意力网络表位预测
- Vaxign-ML: 保护性抗原预测和免疫原性评估
"""

__version__ = "1.0.0"
__author__ = "BacterialAntigenFinder Team"

from .preprocessor import FastaParser, PdbValidator, DataMapper
from .predictors import (
    BasePredictor,
    BepipredWrapper,
    DiscotopeWrapper,
    GraphbepiWrapper,
    EpigraphWrapper,
    VaxignmlWrapper
)
from .aggregator import ConsensusScorer, ImmunogenicityEvaluator, CandidateRanker
from .reporter import CsvExporter, HtmlReporter
