"""
耐药细菌抗原筛选测试

测试Pseudomonas aeruginosa和Klebsiella pneumoniae的抗原筛选功能，
包括FASTA解析、物种配置加载、表位级别排序、CSV导出、VaxiJen抗原性计算
以及端到端Pipeline组件测试。
"""

import os
import sys
import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch
from dataclasses import dataclass

import pytest
import pandas as pd
import numpy as np

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.preprocessor.fasta_parser import FastaParser, ProteinSequence
from src.preprocessor.pdb_validator import PdbValidator
from src.preprocessor.data_mapper import DataMapper
from src.aggregator.consensus_scorer import ConsensusScorer
from src.aggregator.immunogenicity import ImmunogenicityEvaluator
from src.aggregator.candidate_ranker import CandidateRanker
from src.aggregator.vaxijen_calculator import VaxiJenCalculator, VaxiJenResult
from src.reporter.csv_exporter import CsvExporter
from src.reporter.html_report import HtmlReporter
from src.predictors.base_predictor import PredictionResult, EpitopePrediction


# ============================================================
# Fixtures
# ============================================================

PROJECT_ROOT = Path(__file__).parent.parent
EXAMPLE_DATA_DIR = PROJECT_ROOT / "example_data"


@pytest.fixture
def pa_fasta_path():
    """铜绿假单胞菌FASTA文件路径"""
    path = EXAMPLE_DATA_DIR / "pseudomonas_aeruginosa_antigens.fasta"
    if not path.exists():
        pytest.skip(f"测试数据文件不存在: {path}")
    return str(path)


@pytest.fixture
def kp_fasta_path():
    """肺炎克雷伯菌FASTA文件路径"""
    path = EXAMPLE_DATA_DIR / "klebsiella_pneumoniae_antigens.fasta"
    if not path.exists():
        pytest.skip(f"测试数据文件不存在: {path}")
    return str(path)


@pytest.fixture
def pa_sequences(pa_fasta_path):
    """解析铜绿假单胞菌FASTA序列"""
    parser = FastaParser()
    return parser.parse(pa_fasta_path)


@pytest.fixture
def kp_sequences(kp_fasta_path):
    """解析肺炎克雷伯菌FASTA序列"""
    parser = FastaParser()
    return parser.parse(kp_fasta_path)


@pytest.fixture
def default_config():
    """加载默认配置"""
    import yaml
    config_path = PROJECT_ROOT / "config" / "default_config.yaml"
    if not config_path.exists():
        pytest.skip(f"配置文件不存在: {config_path}")
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


@pytest.fixture
def mock_predictions():
    """创建模拟的预测结果，用于测试聚合器和报告器"""
    # 为两个蛋白质创建模拟预测
    proteins = ['OprF_PSEAE', 'OprI_PSEAE']
    predictions = {}

    for model_name, weight in [('bepipred', 0.25), ('discotope', 0.25),
                                ('graphbepi', 0.20), ('epigraph', 0.20),
                                ('vaxignml', 0.10)]:
        preds = []
        for protein_id in proteins:
            for i in range(1, 51):
                # 模拟表位模式：某些区域分数高
                is_epi = (5 <= i <= 15) or (30 <= i <= 40)
                score = np.random.uniform(0.6, 0.95) if is_epi else np.random.uniform(0.1, 0.4)
                preds.append(EpitopePrediction(
                    protein_id=protein_id,
                    residue_id=i,
                    residue_name='ACDEFGHIKLMNPQRSTVWY'[i % 20],
                    score=score,
                    is_epitope=is_epi,
                    additional_info={'subcellular_location': 'OuterMembrane'}
                ))

        predictions[model_name] = PredictionResult(
            predictor_name=model_name,
            predictions=preds,
            protein_scores={pid: np.random.uniform(40, 80) for pid in proteins}
        )

    return predictions


@pytest.fixture
def mock_consensus_df():
    """创建模拟的共识评分DataFrame"""
    data = []
    proteins = ['OprF_PSEAE', 'OprI_PSEAE']

    for protein_id in proteins:
        for i in range(1, 51):
            is_epi = (5 <= i <= 15) or (30 <= i <= 40)
            data.append({
                'protein_id': protein_id,
                'residue_id': i,
                'residue_name': 'ACDEFGHIKLMNPQRSTVWY'[i % 20],
                'consensus_score': np.random.uniform(0.6, 0.9) if is_epi else np.random.uniform(0.1, 0.35),
                'vote_count': 4 if is_epi else 1,
                'is_consensus_epitope': is_epi,
                'bepipred_score': np.random.uniform(0.5, 0.9),
                'bepipred_epitope': is_epi,
                'discotope_score': np.random.uniform(0.5, 0.9),
                'discotope_epitope': is_epi,
                'graphbepi_score': np.random.uniform(0.5, 0.9),
                'graphbepi_epitope': is_epi,
                'epigraph_score': np.random.uniform(0.5, 0.9),
                'epigraph_epitope': is_epi,
                'vaxignml_score': np.random.uniform(0.5, 0.9),
                'vaxignml_epitope': is_epi,
            })

    return pd.DataFrame(data)


@pytest.fixture
def mock_immunogenicity_df():
    """创建模拟的免疫原性评估DataFrame"""
    return pd.DataFrame({
        'protein_id': ['OprF_PSEAE', 'OprI_PSEAE'],
        'protegenicity_score': [75.0, 60.0],
        'is_protective_antigen': [True, True],
        'subcellular_location': ['OuterMembrane', 'OuterMembrane'],
        'has_signal_peptide': [True, True],
        'transmembrane_regions': [1, 0],
        'antigenicity_score': [0.55, 0.48],
        'immunogenicity_rank': [1, 2],
        'recommendation': ['HIGH', 'MEDIUM'],
        'composite_score': [0.72, 0.58]
    })


@pytest.fixture
def sample_sequences():
    """创建用于VaxiJen测试的样本序列"""
    return {
        'test_outer_membrane': ProteinSequence(
            id='test_outer_membrane',
            name='test_outer_membrane',
            description='Test outer membrane protein',
            sequence='MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAG',
            length=60
        ),
        'test_toxin': ProteinSequence(
            id='test_toxin',
            name='test_toxin',
            description='Test exotoxin',
            sequence='MHLTPHWIPLVASLGLLAGGSFASAAEEAFDLWNECAKACVLDLKDGVRSSRMSVDP',
            length=60
        )
    }


# ============================================================
# 测试类: FASTA解析 - 铜绿假单胞菌
# ============================================================

class TestPseudomonasAeruginosaFasta:
    """测试铜绿假单胞菌FASTA文件解析"""

    def test_parse_pa_fasta(self, pa_sequences):
        """测试解析铜绿假单胞菌FASTA文件"""
        assert len(pa_sequences) == 10, f"预期10个序列，实际{len(pa_sequences)}个"

    def test_pa_protein_ids(self, pa_sequences):
        """测试铜绿假单胞菌蛋白质ID"""
        expected_ids = [
            'OprF_PSEAE', 'OprI_PSEAE', 'OprD_PSEAE', 'PilA_PSEAE',
            'PcrV_PSEAE', 'ToxA_PSEAE', 'LasB_PSEAE', 'PhoP_PSEAE',
            'AlgD_PSEAE', 'Azurin_PSEAE'
        ]
        for pid in expected_ids:
            assert pid in pa_sequences, f"缺少蛋白质ID: {pid}"

    def test_pa_sequence_lengths(self, pa_sequences):
        """测试铜绿假单胞菌序列长度"""
        # OprF: 350aa
        assert pa_sequences['OprF_PSEAE'].length == 350
        # OprI: 83aa (最短)
        assert pa_sequences['OprI_PSEAE'].length == 83
        # ToxA: 638aa (最长)
        assert pa_sequences['ToxA_PSEAE'].length == 638

    def test_pa_sequences_valid(self, pa_sequences):
        """测试铜绿假单胞菌序列有效性"""
        for pid, seq in pa_sequences.items():
            assert seq.validate(), f"序列 {pid} 包含非标准氨基酸"

    def test_pa_sequence_content(self, pa_sequences):
        """测试铜绿假单胞菌序列内容"""
        # 验证OprF序列以甲硫氨酸开头
        assert pa_sequences['OprF_PSEAE'].sequence.startswith('M')
        # 验证所有序列都是大写
        for pid, seq in pa_sequences.items():
            assert seq.sequence == seq.sequence.upper(), f"序列 {pid} 不全为大写"

    def test_pa_stats(self, pa_sequences):
        """测试铜绿假单胞菌序列统计"""
        parser = FastaParser()
        stats = parser.get_sequence_stats(pa_sequences)
        assert stats['count'] == 10
        assert stats['total_residues'] > 0
        assert stats['min_length'] == 83
        assert stats['max_length'] == 638


# ============================================================
# 测试类: FASTA解析 - 肺炎克雷伯菌
# ============================================================

class TestKlebsiellaPneumoniaeFasta:
    """测试肺炎克雷伯菌FASTA文件解析"""

    def test_parse_kp_fasta(self, kp_sequences):
        """测试解析肺炎克雷伯菌FASTA文件"""
        assert len(kp_sequences) == 10, f"预期10个序列，实际{len(kp_sequences)}个"

    def test_kp_protein_ids(self, kp_sequences):
        """测试肺炎克雷伯菌蛋白质ID"""
        expected_ids = [
            'OmpK36_KLEPN', 'OmpK35_KLEPN', 'OmpA_KLEPN', 'FimH_KLEPN',
            'MrkD_KLEPN', 'YbtQ_KLEPN', 'EntB_KLEPN', 'KpnO_KLEPN',
            'Peg344_KLEPN', 'RmpA_KLEPN'
        ]
        for pid in expected_ids:
            assert pid in kp_sequences, f"缺少蛋白质ID: {pid}"

    def test_kp_sequence_lengths(self, kp_sequences):
        """测试肺炎克雷伯菌序列长度"""
        # OmpK36: 349aa
        assert kp_sequences['OmpK36_KLEPN'].length == 349
        # RmpA: 217aa
        assert kp_sequences['RmpA_KLEPN'].length == 217
        # YbtQ: 710aa (最长)
        assert kp_sequences['YbtQ_KLEPN'].length == 710

    def test_kp_sequences_valid(self, kp_sequences):
        """测试肺炎克雷伯菌序列有效性"""
        for pid, seq in kp_sequences.items():
            assert seq.validate(), f"序列 {pid} 包含非标准氨基酸"

    def test_kp_stats(self, kp_sequences):
        """测试肺炎克雷伯菌序列统计"""
        parser = FastaParser()
        stats = parser.get_sequence_stats(kp_sequences)
        assert stats['count'] == 10
        assert stats['total_residues'] > 0
        assert stats['max_length'] == 710


# ============================================================
# 测试类: 物种配置加载
# ============================================================

class TestSpeciesConfig:
    """测试物种配置加载"""

    def test_config_has_predefined_species(self, default_config):
        """测试配置文件包含预定义物种"""
        predefined = default_config.get('species', {}).get('predefined', {})
        assert 'Pseudomonas_aeruginosa' in predefined
        assert 'Klebsiella_pneumoniae' in predefined

    def test_pa_species_config(self, default_config):
        """测试铜绿假单胞菌配置"""
        pa_config = default_config['species']['predefined']['Pseudomonas_aeruginosa']
        assert pa_config['organism_type'] == 'gram-'
        assert pa_config['drug_resistance'] == 'MDR/XDR'
        assert 'display_name' in pa_config

    def test_kp_species_config(self, default_config):
        """测试肺炎克雷伯菌配置"""
        kp_config = default_config['species']['predefined']['Klebsiella_pneumoniae']
        assert kp_config['organism_type'] == 'gram-'
        assert kp_config['drug_resistance'] == 'ESBL/CRE'
        assert 'display_name' in kp_config

    def test_all_predefined_are_gram_negative_or_positive(self, default_config):
        """测试所有预定义物种的organism_type有效"""
        predefined = default_config['species']['predefined']
        valid_types = {'gram+', 'gram-'}
        for species_name, config in predefined.items():
            assert config['organism_type'] in valid_types, \
                f"物种 {species_name} 的 organism_type 无效: {config['organism_type']}"

    def test_config_has_model_weights(self, default_config):
        """测试配置文件包含模型权重"""
        weights = default_config.get('models', {}).get('weights', {})
        expected_models = ['bepipred', 'discotope', 'graphbepi', 'epigraph', 'vaxignml']
        for model in expected_models:
            assert model in weights, f"缺少模型权重: {model}"

    def test_model_weights_sum(self, default_config):
        """测试模型权重之和接近1.0"""
        weights = default_config.get('models', {}).get('weights', {})
        total = sum(weights.values())
        assert abs(total - 1.0) < 0.01, f"模型权重之和 {total} 不等于1.0"

    def test_config_has_vaxijen_settings(self, default_config):
        """测试配置文件包含VaxiJen设置"""
        immuno = default_config.get('immunogenicity', {})
        assert 'use_vaxijen' in immuno
        assert immuno['use_vaxijen'] is True
        vaxijen = immuno.get('vaxijen', {})
        assert 'thresholds' in vaxijen
        assert 'bacteria' in vaxijen['thresholds']


# ============================================================
# 测试类: 表位级别排序
# ============================================================

class TestEpitopeRanking:
    """测试表位级别排序功能"""

    def test_rank_epitopes_returns_dataframe(self, mock_consensus_df,
                                              mock_immunogenicity_df,
                                              pa_sequences):
        """测试rank_epitopes返回DataFrame"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        result = ranker.rank_epitopes(
            mock_consensus_df, mock_immunogenicity_df, pa_sequences
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_rank_epitopes_has_required_columns(self, mock_consensus_df,
                                                  mock_immunogenicity_df,
                                                  pa_sequences):
        """测试rank_epitopes输出包含必要列"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        result = ranker.rank_epitopes(
            mock_consensus_df, mock_immunogenicity_df, pa_sequences
        )
        required_cols = [
            'protein_id', 'epitope_start', 'epitope_end',
            'epitope_sequence', 'epitope_length',
            'avg_consensus_score', 'composite_score', 'recommendation', 'rank'
        ]
        for col in required_cols:
            assert col in result.columns, f"缺少列: {col}"

    def test_rank_epitopes_sorted_by_composite(self, mock_consensus_df,
                                                 mock_immunogenicity_df,
                                                 pa_sequences):
        """测试表位按composite_score降序排列"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        result = ranker.rank_epitopes(
            mock_consensus_df, mock_immunogenicity_df, pa_sequences
        )
        if len(result) > 1:
            scores = result['composite_score'].tolist()
            assert scores == sorted(scores, reverse=True), \
                "表位未按composite_score降序排列"

    def test_rank_epitopes_recommendation_valid(self, mock_consensus_df,
                                                  mock_immunogenicity_df,
                                                  pa_sequences):
        """测试推荐等级值有效"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        result = ranker.rank_epitopes(
            mock_consensus_df, mock_immunogenicity_df, pa_sequences
        )
        valid_recs = {'HIGH', 'MEDIUM', 'LOW'}
        for rec in result['recommendation'].unique():
            assert rec in valid_recs, f"无效推荐等级: {rec}"

    def test_rank_epitopes_composite_score_range(self, mock_consensus_df,
                                                   mock_immunogenicity_df,
                                                   pa_sequences):
        """测试composite_score在有效范围内"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        result = ranker.rank_epitopes(
            mock_consensus_df, mock_immunogenicity_df, pa_sequences
        )
        assert (result['composite_score'] >= 0).all()
        assert (result['composite_score'] <= 1).all()

    def test_determine_epitope_recommendation(self):
        """测试表位推荐等级判定逻辑"""
        ranker = CandidateRanker()

        # HIGH: composite >= 0.6, consensus >= 0.5, protegenicity >= 50
        assert ranker._determine_epitope_recommendation(0.7, 0.6, 60, 0.5) == 'HIGH'

        # MEDIUM: composite >= 0.4 and (consensus >= 0.3 or protegenicity >= 30)
        assert ranker._determine_epitope_recommendation(0.5, 0.4, 40, 0.4) == 'MEDIUM'

        # LOW
        assert ranker._determine_epitope_recommendation(0.2, 0.1, 10, 0.2) == 'LOW'

    def test_rank_epitopes_top_n_limit(self, mock_consensus_df,
                                        mock_immunogenicity_df,
                                        pa_sequences):
        """测试top_n限制"""
        ranker = CandidateRanker(top_n=3, min_epitope_length=3)
        result = ranker.rank_epitopes(
            mock_consensus_df, mock_immunogenicity_df, pa_sequences
        )
        assert len(result) <= 3


# ============================================================
# 测试类: CSV导出
# ============================================================

class TestCsvExportWithSpecies:
    """测试带物种信息的CSV导出"""

    def test_export_with_species(self, tmp_path):
        """测试CSV导出包含物种信息"""
        candidates = pd.DataFrame({
            'rank': [1, 2],
            'protein_id': ['OprF_PSEAE', 'OprI_PSEAE'],
            'sequence_length': [350, 83],
            'residue_range': ['5-15', '30-40'],
            'epitope_sequence': ['ACDEFGHIKLMN', 'MNPQRSTVWYAC'],
            'num_epitope_regions': [2, 1],
            'total_epitope_residues': [20, 10],
            'avg_consensus_score': [0.8, 0.7],
            'max_consensus_score': [0.95, 0.85],
            'protegenicity_score': [75.0, 60.0],
            'immunogenicity_rank': [1, 2],
            'subcellular_location': ['OuterMembrane', 'OuterMembrane'],
            'recommendation': ['HIGH', 'MEDIUM'],
            'composite_score': [0.72, 0.58]
        })

        exporter = CsvExporter()
        output_path = tmp_path / 'candidates.csv'
        result_path = exporter.export(candidates, output_path, species='Pseudomonas_aeruginosa')

        assert os.path.exists(result_path)

        # 读取并验证
        with open(result_path, 'r', encoding='utf-8-sig') as f:
            content = f.read()

        # 验证species注释行
        assert '# Species: Pseudomonas_aeruginosa' in content

        # 验证Species列
        df = pd.read_csv(result_path, comment='#')
        assert 'Species' in df.columns
        assert df['Species'].iloc[0] == 'Pseudomonas_aeruginosa'

    def test_export_without_species(self, tmp_path):
        """测试CSV导出不带物种信息"""
        candidates = pd.DataFrame({
            'rank': [1],
            'protein_id': ['test_prot'],
            'sequence_length': [100],
            'residue_range': ['1-10'],
            'epitope_sequence': ['ACDEFGHIKL'],
            'num_epitope_regions': [1],
            'total_epitope_residues': [10],
            'avg_consensus_score': [0.8],
            'max_consensus_score': [0.9],
            'protegenicity_score': [70.0],
            'immunogenicity_rank': [1],
            'subcellular_location': ['OuterMembrane'],
            'recommendation': ['HIGH'],
            'composite_score': [0.75]
        })

        exporter = CsvExporter()
        output_path = tmp_path / 'candidates.csv'
        result_path = exporter.export(candidates, output_path, species=None)

        df = pd.read_csv(result_path)
        assert 'Species' not in df.columns

    def test_export_column_renaming(self, tmp_path):
        """测试CSV导出列名重命名"""
        candidates = pd.DataFrame({
            'rank': [1],
            'protein_id': ['OprF_PSEAE'],
            'sequence_length': [350],
            'avg_consensus_score': [0.8],
            'recommendation': ['HIGH'],
            'composite_score': [0.72]
        })

        exporter = CsvExporter()
        output_path = tmp_path / 'candidates.csv'
        exporter.export(candidates, output_path, species='Pseudomonas_aeruginosa')

        df = pd.read_csv(output_path, comment='#')
        # 验证重命名后的列
        assert 'Protein_ID' in df.columns
        assert 'Avg_Consensus_Score' in df.columns
        assert 'Recommendation' in df.columns


# ============================================================
# 测试类: VaxiJen抗原性计算
# ============================================================

class TestVaxiJenCalculation:
    """测试VaxiJen抗原性计算"""

    def test_vaxijen_calculate_bacteria(self, sample_sequences):
        """测试细菌模型的VaxiJen计算"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        for seq_id, seq in sample_sequences.items():
            result = calc.calculate(seq)
            assert isinstance(result, VaxiJenResult)
            assert 0.0 <= result.antigenicity_score <= 1.0
            assert result.organism_type == 'bacteria'
            assert result.threshold == 0.4

    def test_vaxijen_is_probable_antigen(self, sample_sequences):
        """测试VaxiJen抗原判定"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        for seq_id, seq in sample_sequences.items():
            result = calc.calculate(seq)
            if result.antigenicity_score >= 0.4:
                assert result.is_probable_antigen is True
            else:
                assert result.is_probable_antigen is False

    def test_vaxijen_different_organisms(self, sample_sequences):
        """测试不同生物类型的VaxiJen计算"""
        seq = list(sample_sequences.values())[0]

        calc_bacteria = VaxiJenCalculator(organism_type='bacteria')
        calc_virus = VaxiJenCalculator(organism_type='virus')

        result_b = calc_bacteria.calculate(seq)
        result_v = calc_virus.calculate(seq)

        # 两种模型都应返回有效分数
        assert 0.0 <= result_b.antigenicity_score <= 1.0
        assert 0.0 <= result_v.antigenicity_score <= 1.0

        # 病毒模型阈值也是0.4
        assert result_v.threshold == 0.4

    def test_vaxijen_custom_threshold(self, sample_sequences):
        """测试自定义阈值"""
        calc = VaxiJenCalculator(organism_type='bacteria', custom_threshold=0.6)
        seq = list(sample_sequences.values())[0]
        result = calc.calculate(seq)
        assert result.threshold == 0.6

    def test_vaxijen_batch_calculate(self, sample_sequences):
        """测试批量计算"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        results = calc.batch_calculate(sample_sequences)

        assert len(results) == len(sample_sequences)
        for seq_id, result in results.items():
            assert isinstance(result, VaxiJenResult)
            assert result.protein_id == seq_id

    def test_vaxijen_detailed_analysis(self, sample_sequences):
        """测试详细分析"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        seq = list(sample_sequences.values())[0]
        analysis = calc.get_detailed_analysis(seq)

        assert 'protein_id' in analysis
        assert 'antigenicity_score' in analysis
        assert 'component_scores' in analysis
        assert 'aa_composition' in analysis
        assert 'interpretation' in analysis

        # 验证各分量分数
        components = analysis['component_scores']
        for comp_name in ['acc_score', 'composition_score', 'hydrophobicity_score',
                          'hydrophilicity_score', 'diversity_score']:
            assert comp_name in components
            assert 0.0 <= components[comp_name] <= 1.0

    def test_vaxijen_pa_sequences(self, pa_sequences):
        """测试VaxiJen对铜绿假单胞菌序列的计算"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        results = calc.batch_calculate(pa_sequences)

        assert len(results) == 10
        # 外膜蛋白应该有较高的抗原性
        oprf_result = results.get('OprF_PSEAE')
        if oprf_result:
            assert oprf_result.antigenicity_score > 0, "OprF抗原性分数应为正值"

    def test_vaxijen_kp_sequences(self, kp_sequences):
        """测试VaxiJen对肺炎克雷伯菌序列的计算"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        results = calc.batch_calculate(kp_sequences)

        assert len(results) == 10
        # 所有结果应有效
        for seq_id, result in results.items():
            assert 0.0 <= result.antigenicity_score <= 1.0

    def test_vaxijen_short_sequence(self):
        """测试VaxiJen对短序列的处理"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        short_seq = ProteinSequence(
            id='short', name='short', description='',
            sequence='ACDE', length=4
        )
        result = calc.calculate(short_seq)
        assert result.antigenicity_score == 0.0
        assert result.is_probable_antigen is False

    def test_vaxijen_acc_descriptors(self, sample_sequences):
        """测试ACC描述符计算"""
        calc = VaxiJenCalculator(organism_type='bacteria')
        seq = list(sample_sequences.values())[0]
        result = calc.calculate(seq)

        if result.acc_descriptors is not None:
            # ACC描述符应为 5 * lag 维度
            expected_dim = 5 * calc.lag
            assert len(result.acc_descriptors) == expected_dim


# ============================================================
# 测试类: 端到端Pipeline组件（Mock预测器）
# ============================================================

class TestPipelineComponents:
    """测试Pipeline各组件的集成"""

    def test_consensus_scorer_with_mock_predictions(self, mock_predictions):
        """测试共识评分器与模拟预测结果"""
        scorer = ConsensusScorer(
            weights={'bepipred': 0.25, 'discotope': 0.25,
                     'graphbepi': 0.20, 'epigraph': 0.20, 'vaxignml': 0.10},
            threshold=0.5,
            min_votes=2
        )

        result = scorer.compute_consensus(mock_predictions)

        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert 'consensus_score' in result.columns
        assert 'is_consensus_epitope' in result.columns
        assert 'vote_count' in result.columns

    def test_candidate_ranker_with_mock_data(self, mock_consensus_df,
                                              mock_immunogenicity_df,
                                              pa_sequences):
        """测试候选排序器与模拟数据"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        result = ranker.rank(mock_consensus_df, mock_immunogenicity_df, pa_sequences)

        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert 'rank' in result.columns
        assert 'protein_id' in result.columns
        assert 'composite_score' in result.columns
        assert 'recommendation' in result.columns

    def test_csv_exporter_with_ranked_candidates(self, mock_consensus_df,
                                                   mock_immunogenicity_df,
                                                   pa_sequences, tmp_path):
        """测试CSV导出器与排序后的候选"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        candidates = ranker.rank(mock_consensus_df, mock_immunogenicity_df, pa_sequences)

        exporter = CsvExporter()
        output_path = tmp_path / 'antigen_candidates.csv'
        result_path = exporter.export(candidates, output_path, species='Pseudomonas_aeruginosa')

        assert os.path.exists(result_path)
        df = pd.read_csv(result_path, comment='#')
        assert len(df) > 0
        assert 'Protein_ID' in df.columns

    def test_html_reporter_with_mock_results(self, mock_consensus_df,
                                               mock_immunogenicity_df,
                                               pa_sequences, tmp_path):
        """测试HTML报告生成器"""
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        candidates = ranker.rank(mock_consensus_df, mock_immunogenicity_df, pa_sequences)

        results = {
            'candidates': candidates,
            'consensus': mock_consensus_df,
            'predictions': {},
            'immunogenicity': mock_immunogenicity_df
        }

        reporter = HtmlReporter()
        output_path = tmp_path / 'analysis_report.html'
        result_path = reporter.generate(results, output_path, species='Pseudomonas_aeruginosa')

        assert os.path.exists(result_path)
        with open(result_path, 'r', encoding='utf-8') as f:
            content = f.read()
        assert 'BacterialAntigenFinder' in content
        assert 'Pseudomonas_aeruginosa' in content

    def test_full_mock_pipeline(self, pa_sequences, tmp_path):
        """测试完整的模拟Pipeline流程"""
        # Step 1: 序列已解析 (fixture提供)

        # Step 2: 创建模拟预测
        mock_preds = {}
        for model_name in ['bepipred', 'discotope', 'graphbepi', 'epigraph', 'vaxignml']:
            preds = []
            for protein_id in list(pa_sequences.keys())[:3]:
                seq = pa_sequences[protein_id]
                for i in range(1, min(seq.length + 1, 51)):
                    is_epi = (5 <= i <= 15) or (30 <= i <= 40)
                    score = 0.8 if is_epi else 0.2
                    preds.append(EpitopePrediction(
                        protein_id=protein_id,
                        residue_id=i,
                        residue_name=seq.sequence[i - 1] if i <= len(seq.sequence) else 'A',
                        score=score,
                        is_epitope=is_epi,
                        additional_info={'subcellular_location': 'OuterMembrane'}
                    ))
            mock_preds[model_name] = PredictionResult(
                predictor_name=model_name,
                predictions=preds,
                protein_scores={pid: 65.0 for pid in list(pa_sequences.keys())[:3]}
            )

        # Step 3: 共识评分
        scorer = ConsensusScorer(threshold=0.5, min_votes=2)
        consensus_df = scorer.compute_consensus(mock_preds)
        assert not consensus_df.empty

        # Step 4: 免疫原性评估（使用mock）
        immuno_df = pd.DataFrame({
            'protein_id': list(pa_sequences.keys())[:3],
            'protegenicity_score': [75.0, 60.0, 55.0],
            'is_protective_antigen': [True, True, True],
            'subcellular_location': ['OuterMembrane', 'OuterMembrane', 'Periplasmic'],
            'has_signal_peptide': [True, True, False],
            'transmembrane_regions': [1, 0, 2],
            'antigenicity_score': [0.55, 0.48, 0.42],
            'immunogenicity_rank': [1, 2, 3],
            'recommendation': ['HIGH', 'MEDIUM', 'MEDIUM'],
            'composite_score': [0.72, 0.58, 0.50]
        })

        # Step 5: 候选排序
        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        candidates = ranker.rank(consensus_df, immuno_df, pa_sequences)
        assert not candidates.empty

        # Step 6: 导出结果
        csv_exporter = CsvExporter()
        csv_path = tmp_path / 'antigen_candidates.csv'
        csv_exporter.export(candidates, csv_path, species='Pseudomonas_aeruginosa')
        assert csv_path.exists()

        # 表位级别排序
        epitope_df = ranker.rank_epitopes(consensus_df, immuno_df, pa_sequences)
        if not epitope_df.empty:
            epitope_csv = tmp_path / 'epitope_candidates.csv'
            epitope_df.to_csv(epitope_csv, index=False)
            assert epitope_csv.exists()

        # HTML报告
        html_reporter = HtmlReporter()
        html_path = tmp_path / 'analysis_report.html'
        html_reporter.generate(
            {'candidates': candidates, 'consensus': consensus_df,
             'predictions': mock_preds, 'immunogenicity': immuno_df},
            html_path,
            species='Pseudomonas_aeruginosa'
        )
        assert html_path.exists()

        # 元数据
        metadata = {
            'version': '1.0.0',
            'species': 'Pseudomonas_aeruginosa',
            'input': {'num_sequences': len(pa_sequences)},
            'output': {'num_candidates': len(candidates)}
        }
        metadata_path = tmp_path / 'run_metadata.json'
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        assert metadata_path.exists()

    def test_immunogenicity_evaluator_with_vaxijen(self, sample_sequences):
        """测试免疫原性评估器使用VaxiJen"""
        evaluator = ImmunogenicityEvaluator(
            protegenicity_threshold=50.0,
            antigenicity_threshold=0.4,
            organism_type='bacteria',
            use_vaxijen=True
        )

        # 创建模拟的Vaxign-ML结果
        mock_vaxignml = PredictionResult(
            predictor_name='vaxignml',
            predictions=[
                EpitopePrediction(
                    protein_id=seq_id, residue_id=1, residue_name='M',
                    score=0.7, is_epitope=True,
                    additional_info={
                        'subcellular_location': 'OuterMembrane',
                        'has_signal_peptide': True,
                        'transmembrane_regions': 1,
                        'protegenicity': 75.0
                    }
                )
                for seq_id in sample_sequences
            ],
            protein_scores={seq_id: 75.0 for seq_id in sample_sequences}
        )

        result = evaluator.evaluate(mock_vaxignml, sample_sequences)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == len(sample_sequences)
        assert 'antigenicity_score' in result.columns
        assert 'recommendation' in result.columns

        # VaxiJen计算的抗原性分数应大于0
        assert (result['antigenicity_score'] > 0).all()


# ============================================================
# 测试类: 跨物种对比
# ============================================================

class TestCrossSpeciesComparison:
    """测试跨物种对比功能"""

    def test_both_species_parseable(self, pa_sequences, kp_sequences):
        """测试两种菌的FASTA文件都能正确解析"""
        assert len(pa_sequences) == 10
        assert len(kp_sequences) == 10

    def test_species_protein_ids_dont_overlap(self, pa_sequences, kp_sequences):
        """测试两种菌的蛋白质ID不重叠"""
        pa_ids = set(pa_sequences.keys())
        kp_ids = set(kp_sequences.keys())
        overlap = pa_ids & kp_ids
        assert len(overlap) == 0, f"蛋白质ID重叠: {overlap}"

    def test_vaxijen_both_species(self, pa_sequences, kp_sequences):
        """测试VaxiJen对两种菌的抗原性计算"""
        calc = VaxiJenCalculator(organism_type='bacteria')

        pa_results = calc.batch_calculate(pa_sequences)
        kp_results = calc.batch_calculate(kp_sequences)

        # 两种菌都应有10个结果
        assert len(pa_results) == 10
        assert len(kp_results) == 10

        # 所有分数应在有效范围内
        for results in [pa_results, kp_results]:
            for result in results.values():
                assert 0.0 <= result.antigenicity_score <= 1.0

    def test_ranker_with_both_species(self, mock_consensus_df,
                                       mock_immunogenicity_df,
                                       pa_sequences, kp_sequences, tmp_path):
        """测试排序器对两种菌的序列"""
        # 合并序列
        all_sequences = {**pa_sequences, **kp_sequences}

        ranker = CandidateRanker(top_n=50, min_epitope_length=3)
        candidates = ranker.rank(mock_consensus_df, mock_immunogenicity_df, all_sequences)

        # 导出带物种信息
        exporter = CsvExporter()
        csv_path = tmp_path / 'combined_candidates.csv'
        exporter.export(candidates, csv_path, species='Pseudomonas_aeruginosa')
        assert csv_path.exists()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
