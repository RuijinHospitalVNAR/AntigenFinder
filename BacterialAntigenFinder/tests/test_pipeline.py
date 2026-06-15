"""
Pipeline集成测试

测试整个抗原筛选Pipeline的功能。
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path
import pytest

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.preprocessor.fasta_parser import FastaParser, ProteinSequence
from src.preprocessor.pdb_validator import PdbValidator, ProteinStructure
from src.preprocessor.data_mapper import DataMapper, MappedData
from src.aggregator.consensus_scorer import ConsensusScorer
from src.aggregator.immunogenicity import ImmunogenicityEvaluator
from src.aggregator.candidate_ranker import CandidateRanker
from src.reporter.csv_exporter import CsvExporter
from src.predictors.base_predictor import PredictionResult, EpitopePrediction


class TestFastaParser:
    """测试FASTA解析器"""
    
    def test_parse_valid_fasta(self, tmp_path):
        """测试解析有效的FASTA文件"""
        # 创建测试文件
        fasta_content = """>protein1 Test protein 1
MKFLILLFNILCLFPVLAADNHGVGPQGASGVDPITFDINSNQTGVQSLTFDA
>protein2 Test protein 2
MKKLLILTLLFGIAGPAIAAQYEEVVNNNGPTHENQLGAG
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)
        
        parser = FastaParser()
        sequences = parser.parse(str(fasta_file))
        
        assert len(sequences) == 2
        assert 'protein1' in sequences
        assert 'protein2' in sequences
        assert sequences['protein1'].length == 53
        assert sequences['protein2'].length == 40
    
    def test_validate_sequence(self):
        """测试序列验证"""
        valid_seq = ProteinSequence(
            id='test',
            name='test',
            description='test protein',
            sequence='ACDEFGHIKLMNPQRSTVWY',
            length=20
        )
        assert valid_seq.validate() == True
        
        # 含非标准氨基酸的序列仍然有效（跳过非标准字符）
        seq_with_x = ProteinSequence(
            id='test',
            name='test',
            description='test protein',
            sequence='ACDEFXGHIK',
            length=10
        )
        # X不在标准氨基酸中，validate会返回False
        assert seq_with_x.validate() == False
    
    def test_filter_by_length(self):
        """测试按长度过滤"""
        sequences = {
            'short': ProteinSequence('short', 'short', '', 'ACDEF', 5),
            'medium': ProteinSequence('medium', 'medium', '', 'A' * 100, 100),
            'long': ProteinSequence('long', 'long', '', 'A' * 6000, 6000)
        }
        
        parser = FastaParser()
        filtered = parser.filter_by_length(sequences, min_length=50, max_length=5000)
        
        assert len(filtered) == 1
        assert 'medium' in filtered


class TestConsensusScorer:
    """测试共识评分器"""
    
    def test_weighted_average(self):
        """测试加权平均计算"""
        scorer = ConsensusScorer(
            weights={'model1': 0.6, 'model2': 0.4},
            threshold=0.5
        )
        
        scores = {'model1': 0.8, 'model2': 0.6}
        result = scorer._weighted_average(scores)
        
        expected = (0.8 * 0.6 + 0.6 * 0.4) / 1.0
        assert abs(result - expected) < 0.001
    
    def test_compute_consensus(self):
        """测试共识计算"""
        # 创建模拟预测结果
        predictions = {
            'model1': PredictionResult(
                predictor_name='model1',
                predictions=[
                    EpitopePrediction('prot1', 1, 'A', 0.8, True),
                    EpitopePrediction('prot1', 2, 'C', 0.3, False),
                ],
                protein_scores={'prot1': 0.55}
            ),
            'model2': PredictionResult(
                predictor_name='model2',
                predictions=[
                    EpitopePrediction('prot1', 1, 'A', 0.7, True),
                    EpitopePrediction('prot1', 2, 'C', 0.4, False),
                ],
                protein_scores={'prot1': 0.55}
            )
        }
        
        scorer = ConsensusScorer(
            weights={'model1': 0.5, 'model2': 0.5},
            threshold=0.5,
            min_votes=2
        )
        
        result = scorer.compute_consensus(predictions)
        
        assert len(result) == 2
        assert result[result['residue_id'] == 1]['is_consensus_epitope'].values[0] == True
        assert result[result['residue_id'] == 2]['is_consensus_epitope'].values[0] == False


class TestCandidateRanker:
    """测试候选排序器"""
    
    def test_extract_epitope_regions(self):
        """测试表位区域提取"""
        import pandas as pd
        
        # 创建模拟数据
        data = {
            'protein_id': ['prot1'] * 10,
            'residue_id': list(range(1, 11)),
            'residue_name': list('ACDEFGHIKL'),
            'consensus_score': [0.8, 0.7, 0.6, 0.65, 0.7, 0.3, 0.2, 0.8, 0.75, 0.8],
            'is_consensus_epitope': [True, True, True, True, True, False, False, True, True, True]
        }
        df = pd.DataFrame(data)
        
        ranker = CandidateRanker(min_epitope_length=3)
        regions = ranker._extract_epitope_regions(df)
        
        # 应该有两个区域：1-5 和 8-10
        assert len(regions) == 2
        assert regions[0][0] == 1  # 第一个区域起始
        assert regions[0][1] == 5  # 第一个区域结束


class TestDataMapper:
    """测试数据映射器"""
    
    def test_normalize_id(self):
        """测试ID规范化"""
        mapper = DataMapper()
        
        assert mapper._normalize_id('protein_A.pdb') == 'protein'
        assert mapper._normalize_id('PROTEIN_alphafold') == 'protein'
        assert mapper._normalize_id('test_1') == 'test'
    
    def test_calculate_sequence_similarity(self):
        """测试序列相似度计算"""
        mapper = DataMapper()
        
        similarity = mapper._calculate_sequence_similarity('ACDEF', 'ACDEF')
        assert similarity == 1.0
        
        similarity = mapper._calculate_sequence_similarity('ACDEF', 'ACXYZ')
        assert similarity == 0.4  # 2/5


class TestCsvExporter:
    """测试CSV导出器"""
    
    def test_export(self, tmp_path):
        """测试CSV导出"""
        import pandas as pd
        
        candidates = pd.DataFrame({
            'rank': [1, 2],
            'protein_id': ['prot1', 'prot2'],
            'residue_range': ['1-10', '20-30'],
            'epitope_sequence': ['ACDEFGHIKL', 'MNPQRSTVWY'],
            'avg_consensus_score': [0.8, 0.7],
            'recommendation': ['HIGH', 'MEDIUM']
        })
        
        exporter = CsvExporter()
        output_path = tmp_path / 'candidates.csv'
        result_path = exporter.export(candidates, output_path)
        
        assert os.path.exists(result_path)
        
        # 验证内容
        loaded = pd.read_csv(result_path)
        assert len(loaded) == 2


class TestImmunogenicityEvaluator:
    """测试免疫原性评估器"""
    
    def test_calculate_antigenicity(self):
        """测试抗原性计算"""
        evaluator = ImmunogenicityEvaluator()
        
        seq = ProteinSequence(
            id='test',
            name='test',
            description='',
            sequence='AAAAAVVVVV',  # A和V有不同的抗原性系数
            length=10
        )
        
        score = evaluator._calculate_antigenicity(seq)
        
        # 应该是 (5*1.064 + 5*1.383) / 10 = 1.2235
        expected = (5 * 1.064 + 5 * 1.383) / 10
        assert abs(score - expected) < 0.01
    
    def test_get_location_weight(self):
        """测试亚细胞定位权重"""
        evaluator = ImmunogenicityEvaluator()
        
        assert evaluator._get_location_weight('OuterMembrane') == 1.5
        assert evaluator._get_location_weight('Cytoplasmic') == 0.8
        assert evaluator._get_location_weight('Unknown') == 1.0


# 集成测试
class TestIntegration:
    """集成测试"""
    
    @pytest.fixture
    def sample_data(self, tmp_path):
        """创建示例数据"""
        # 创建FASTA文件
        fasta_content = """>OmpA_ECOLI Outer membrane protein A
MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINN
>FimH_ECOLI Type 1 fimbrial adhesin
MKKLLILTLLFGIAGPAIAAQYEEVVNNNGPTHENQLGAG
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)
        
        # 创建PDB目录（空的，仅用于测试）
        pdb_dir = tmp_path / "structures"
        pdb_dir.mkdir()
        
        return {
            'fasta_file': str(fasta_file),
            'pdb_dir': str(pdb_dir),
            'output_dir': str(tmp_path / 'output')
        }
    
    def test_preprocessor_integration(self, sample_data):
        """测试预处理器集成"""
        parser = FastaParser()
        sequences = parser.parse(sample_data['fasta_file'])
        
        assert len(sequences) == 2
        
        # 统计信息
        stats = parser.get_sequence_stats(sequences)
        assert stats['count'] == 2
        assert stats['total_residues'] > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
