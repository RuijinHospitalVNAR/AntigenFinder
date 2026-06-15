"""
CSV导出模块

导出抗原候选清单和详细结果。
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
from datetime import datetime
import pandas as pd


class CsvExporter:
    """CSV导出器"""
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def export(self, 
               candidates_df: pd.DataFrame,
               output_path: Union[str, Path],
               include_metadata: bool = True,
               species: Optional[str] = None) -> str:
        """
        导出抗原候选清单到CSV
        
        Args:
            candidates_df: 候选抗原DataFrame
            output_path: 输出文件路径
            include_metadata: 是否包含元数据
            species: 目标细菌种类名称
        
        Returns:
            输出文件路径
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 准备导出数据
        export_df = candidates_df.copy()
        
        # 添加species列
        if species:
            export_df.insert(0, 'Species', species)
        
        # 重命名列为更易读的名称
        column_mapping = {
            'rank': 'Rank',
            'protein_id': 'Protein_ID',
            'sequence_length': 'Sequence_Length',
            'residue_range': 'Best_Epitope_Range',
            'epitope_sequence': 'Best_Epitope_Sequence',
            'num_epitope_regions': 'Num_Epitope_Regions',
            'total_epitope_residues': 'Total_Epitope_Residues',
            'avg_consensus_score': 'Avg_Consensus_Score',
            'max_consensus_score': 'Max_Consensus_Score',
            'protegenicity_score': 'Protegenicity_Score',
            'immunogenicity_rank': 'Immunogenicity_Rank',
            'subcellular_location': 'Subcellular_Location',
            'recommendation': 'Recommendation',
            'composite_score': 'Composite_Score',
            'all_epitope_regions': 'All_Epitope_Regions'
        }
        
        export_df = export_df.rename(columns=column_mapping)
        
        # 保存CSV
        if species:
            # 写入species信息作为注释头行
            with open(output_path, 'w', encoding='utf-8-sig') as f:
                f.write(f"# Species: {species}\n")
            export_df.to_csv(output_path, mode='a', index=False, encoding='utf-8-sig')
        else:
            export_df.to_csv(output_path, index=False, encoding='utf-8-sig')
        
        self.logger.info(f"已导出 {len(export_df)} 个候选到 {output_path}")
        
        return str(output_path)
    
    def export_detailed_results(self,
                                 consensus_df: pd.DataFrame,
                                 output_dir: Union[str, Path],
                                 prefix: str = 'detailed') -> Dict[str, str]:
        """
        导出详细的残基级别结果
        
        Args:
            consensus_df: 共识结果DataFrame
            output_dir: 输出目录
            prefix: 文件名前缀
        
        Returns:
            文件名到路径的字典
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        exported_files = {}
        
        # 导出完整的残基级别结果
        full_path = output_dir / f'{prefix}_residue_predictions.csv'
        consensus_df.to_csv(full_path, index=False, encoding='utf-8-sig')
        exported_files['full_results'] = str(full_path)
        
        # 导出表位残基
        epitope_df = consensus_df[consensus_df['is_consensus_epitope'] == True]
        epitope_path = output_dir / f'{prefix}_epitope_residues.csv'
        epitope_df.to_csv(epitope_path, index=False, encoding='utf-8-sig')
        exported_files['epitope_residues'] = str(epitope_path)
        
        # 按蛋白质分组导出
        protein_dir = output_dir / 'per_protein'
        protein_dir.mkdir(exist_ok=True)
        
        for protein_id in consensus_df['protein_id'].unique():
            protein_df = consensus_df[consensus_df['protein_id'] == protein_id]
            protein_path = protein_dir / f'{protein_id}_predictions.csv'
            protein_df.to_csv(protein_path, index=False, encoding='utf-8-sig')
        
        exported_files['per_protein_dir'] = str(protein_dir)
        
        self.logger.info(f"已导出详细结果到 {output_dir}")
        
        return exported_files
    
    def export_summary_statistics(self,
                                   results: dict,
                                   output_path: Union[str, Path]) -> str:
        """
        导出统计摘要
        
        Args:
            results: 结果字典
            output_path: 输出文件路径
        
        Returns:
            输出文件路径
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        summary_data = []
        
        # 总体统计
        candidates = results.get('candidates', pd.DataFrame())
        consensus = results.get('consensus', pd.DataFrame())
        
        summary_data.append({
            'Category': 'Overview',
            'Metric': 'Total Proteins Analyzed',
            'Value': consensus['protein_id'].nunique() if not consensus.empty else 0
        })
        
        summary_data.append({
            'Category': 'Overview',
            'Metric': 'Total Residues',
            'Value': len(consensus) if not consensus.empty else 0
        })
        
        summary_data.append({
            'Category': 'Overview',
            'Metric': 'Epitope Residues',
            'Value': consensus['is_consensus_epitope'].sum() if not consensus.empty else 0
        })
        
        summary_data.append({
            'Category': 'Overview',
            'Metric': 'Top Candidates',
            'Value': len(candidates) if not candidates.empty else 0
        })
        
        # 推荐统计
        if not candidates.empty and 'recommendation' in candidates.columns:
            for rec in ['HIGH', 'MEDIUM', 'LOW']:
                count = (candidates['recommendation'] == rec).sum()
                summary_data.append({
                    'Category': 'Recommendations',
                    'Metric': f'{rec} Priority Candidates',
                    'Value': count
                })
        
        # 预测器统计
        predictions = results.get('predictions', {})
        for predictor_name, pred_result in predictions.items():
            if hasattr(pred_result, 'predictions'):
                epitope_count = sum(1 for p in pred_result.predictions if p.is_epitope)
                summary_data.append({
                    'Category': 'Predictors',
                    'Metric': f'{predictor_name} Epitope Predictions',
                    'Value': epitope_count
                })
        
        # 保存
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(output_path, index=False, encoding='utf-8-sig')
        
        self.logger.info(f"已导出统计摘要到 {output_path}")
        
        return str(output_path)
    
    def export_for_visualization(self,
                                  consensus_df: pd.DataFrame,
                                  output_path: Union[str, Path]) -> str:
        """
        导出用于可视化的数据
        
        Args:
            consensus_df: 共识结果DataFrame
            output_path: 输出文件路径
        
        Returns:
            输出文件路径
        """
        output_path = Path(output_path)
        
        # 准备可视化数据
        viz_data = []
        
        for protein_id in consensus_df['protein_id'].unique():
            protein_df = consensus_df[consensus_df['protein_id'] == protein_id]
            
            for _, row in protein_df.iterrows():
                viz_data.append({
                    'protein_id': protein_id,
                    'position': row['residue_id'],
                    'residue': row['residue_name'],
                    'score': row['consensus_score'],
                    'is_epitope': row['is_consensus_epitope'],
                    'votes': row['vote_count']
                })
        
        viz_df = pd.DataFrame(viz_data)
        viz_df.to_csv(output_path, index=False, encoding='utf-8-sig')
        
        return str(output_path)
