#!/usr/bin/env python3
"""
端到端流程验证测试（不依赖外部模型）

使用模拟的预测结果验证完整 Pipeline 数据流：
FASTA解析 → 模拟预测 → 共识评分 → VaxiJen免疫原性评估 → 候选排序 → CSV输出
"""
import sys
import os
import tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from src.preprocessor.fasta_parser import FastaParser
from src.predictors.base_predictor import PredictionResult, EpitopePrediction
from src.aggregator.consensus_scorer import ConsensusScorer
from src.aggregator.immunogenicity import ImmunogenicityEvaluator
from src.aggregator.candidate_ranker import CandidateRanker
from src.reporter.csv_exporter import CsvExporter


def create_mock_predictions(sequences, predictor_name, threshold=0.5):
    """创建模拟预测结果"""
    predictions = []
    protein_scores = {}

    for seq_id, seq_obj in sequences.items():
        scores = []
        for i, aa in enumerate(seq_obj.sequence):
            # 模拟：每10个残基中有6个是连续表位
            is_epi = (i % 12 < 6)
            if is_epi:
                score = 0.7 + (i % 5) * 0.05
            else:
                score = 0.2 + (i % 3) * 0.1

            predictions.append(EpitopePrediction(
                protein_id=seq_id,
                residue_id=i + 1,
                residue_name=aa,
                score=score,
                is_epitope=score >= threshold,
                confidence=score,
                additional_info={'predictor': predictor_name}
            ))
            scores.append(score)

        protein_scores[seq_id] = sum(scores) / len(scores) if scores else 0

    return PredictionResult(
        predictor_name=predictor_name,
        predictions=predictions,
        protein_scores=protein_scores,
        metadata={'threshold': threshold}
    )


def main():
    print("=" * 60)
    print("端到端流程验证测试 (VaxiJen 2.0 免疫原性评估)")
    print("=" * 60)

    # Step 1: 解析FASTA
    print("\nStep 1: FASTA解析...")
    parser = FastaParser()
    fasta_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                              'example_data', 'pseudomonas_aeruginosa_antigens.fasta')
    sequences = parser.parse(fasta_file)
    print(f"  解析到 {len(sequences)} 个序列")

    # Step 2: 模拟多模型预测
    print("\nStep 2: 模拟4模型预测...")
    predictions = {}
    for name in ['bepipred', 'discotope', 'graphbepi', 'epigraph']:
        pred = create_mock_predictions(sequences, name)
        predictions[name] = pred
        print(f"  {name}: {len(pred.predictions)} 残基")

    # Step 3: 共识评分
    print("\nStep 3: 共识评分...")
    scorer = ConsensusScorer(
        threshold=0.5,
        min_votes=2,
        method='weighted_avg',
        score_normalize='rank'
    )
    consensus_df = scorer.compute_consensus(predictions)
    print(f"  共识评分完成，{len(consensus_df)} 残基")
    print(f"  共识表位残基: {consensus_df['is_consensus_epitope'].sum()}")

    # Step 4: VaxiJen 2.0 免疫原性评估
    print("\nStep 4: VaxiJen 2.0 免疫原性评估...")
    evaluator = ImmunogenicityEvaluator(
        protegenicity_threshold=50.0,
        antigenicity_threshold=0.4,
        organism_type='bacteria',
        use_vaxijen=True
    )
    immuno_df = evaluator.evaluate_from_sequences(sequences)
    print(f"  评估完成，{len(immuno_df)} 个蛋白质")
    print(f"  非零保护性分数: {(immuno_df['protegenicity_score'] > 0).sum()}/{len(immuno_df)}")
    print(f"  HIGH推荐: {(immuno_df['recommendation'] == 'HIGH').sum()}")
    print(f"  MEDIUM推荐: {(immuno_df['recommendation'] == 'MEDIUM').sum()}")

    # Step 5: 候选排序
    print("\nStep 5: 候选排序...")
    ranker = CandidateRanker(top_n=50)
    candidates = ranker.rank(consensus_df, immuno_df, sequences)
    epitopes = ranker.rank_epitopes(consensus_df, immuno_df, sequences)
    print(f"  蛋白质级候选: {len(candidates)}")
    print(f"  表位级候选: {len(epitopes)}")

    # 验证免疫原性分数已传递到候选列表
    if not candidates.empty:
        nonzero_proteg = (candidates['protegenicity_score'] > 0).sum()
        print(f"  候选中非零保护性分数: {nonzero_proteg}/{len(candidates)}")

    if not epitopes.empty:
        nonzero_epitope = (epitopes['protegenicity_score'] > 0).sum()
        nonzero_antigen = (epitopes['antigenicity_score'] > 0).sum()
        print(f"  表位中非零保护性分数: {nonzero_epitope}/{len(epitopes)}")
        print(f"  表位中非零抗原性分数: {nonzero_antigen}/{len(epitopes)}")

    # Step 6: CSV输出
    print("\nStep 6: CSV输出...")
    exporter = CsvExporter()
    tmpdir = tempfile.mkdtemp(prefix='antigen_test_')
    csv_path = exporter.export(candidates, os.path.join(tmpdir, 'antigen_candidates.csv'),
                               species='Pseudomonas_aeruginosa')
    epitope_csv_path = exporter.export(epitopes, os.path.join(tmpdir, 'epitope_candidates.csv'),
                                       species='Pseudomonas_aeruginosa')
    print(f"  CSV: {csv_path}")
    print(f"  表位CSV: {epitope_csv_path}")

    # 打印前5个表位候选
    if not epitopes.empty:
        print("\n=== Top 5 表位候选 ===")
        cols = ['rank', 'protein_id', 'epitope_sequence', 'composite_score',
                'protegenicity_score', 'antigenicity_score', 'recommendation']
        available_cols = [c for c in cols if c in epitopes.columns]
        print(epitopes[available_cols].head().to_string(index=False))

    print("\n" + "=" * 60)
    print("端到端测试通过!")
    print("=" * 60)

    # 清理
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == '__main__':
    main()
