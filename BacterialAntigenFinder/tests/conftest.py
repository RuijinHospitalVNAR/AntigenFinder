"""
Pytest配置和共享fixtures
"""

import os
import sys
from pathlib import Path
import tempfile
import pytest

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture
def project_root():
    """返回项目根目录"""
    return Path(__file__).parent.parent


@pytest.fixture
def sample_fasta(tmp_path):
    """创建示例FASTA文件"""
    content = """>protein1 Test protein 1
MKFLILLFNILCLFPVLAADNHGVGPQGASGVDPITFDINSNQTGVQSLTFDA
>protein2 Test protein 2
MKKLLILTLLFGIAGPAIAAQYEEVVNNNGPTHENQLGAG
>protein3 Test protein 3 with longer sequence
MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGA
FGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDVYTRLGG
MVWRADTKSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDN
"""
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(content)
    return str(fasta_file)


@pytest.fixture
def sample_pdb_dir(tmp_path):
    """创建示例PDB目录"""
    pdb_dir = tmp_path / "structures"
    pdb_dir.mkdir()
    
    # 创建一个最小的PDB文件
    pdb_content = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2       3.323   1.513   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2       3.970   2.822   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2       5.481   2.716   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2       6.088   1.643   0.000  1.00 20.00           O
END
"""
    pdb_file = pdb_dir / "protein1.pdb"
    pdb_file.write_text(pdb_content)
    
    return str(pdb_dir)


@pytest.fixture
def output_dir(tmp_path):
    """创建输出目录"""
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    return str(out_dir)


@pytest.fixture
def mock_prediction_result():
    """创建模拟预测结果"""
    from src.predictors.base_predictor import PredictionResult, EpitopePrediction
    
    predictions = [
        EpitopePrediction('protein1', i, 'A', 0.8 if i % 3 == 0 else 0.3, i % 3 == 0)
        for i in range(1, 51)
    ]
    
    return PredictionResult(
        predictor_name='test_predictor',
        predictions=predictions,
        protein_scores={'protein1': 0.5}
    )


@pytest.fixture
def mock_sequences():
    """创建模拟序列字典"""
    from src.preprocessor.fasta_parser import ProteinSequence
    
    return {
        'protein1': ProteinSequence(
            id='protein1',
            name='protein1',
            description='Test protein 1',
            sequence='A' * 50,
            length=50
        ),
        'protein2': ProteinSequence(
            id='protein2',
            name='protein2',
            description='Test protein 2',
            sequence='G' * 40,
            length=40
        )
    }
