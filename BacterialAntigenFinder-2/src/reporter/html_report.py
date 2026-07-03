"""
HTML报告生成模块

生成交互式HTML分析报告。
"""

import logging
import os
from pathlib import Path
from typing import Dict, Optional, Union
from datetime import datetime
import json
import pandas as pd

# 尝试导入plotly，如果不可用则使用简单HTML
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


class HtmlReporter:
    """HTML报告生成器"""
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def generate(self, 
                 results: dict,
                 output_path: Union[str, Path],
                 species: Optional[str] = None) -> str:
        """
        生成HTML分析报告
        
        Args:
            results: 结果字典，包含candidates, consensus, predictions等
            output_path: 输出文件路径
            species: 目标细菌种类名称
        
        Returns:
            输出文件路径
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 获取结果数据
        candidates = results.get('candidates', pd.DataFrame())
        epitope_candidates = results.get('epitope_candidates', pd.DataFrame())
        consensus = results.get('consensus', pd.DataFrame())
        predictions = results.get('predictions', {})
        immunogenicity = results.get('immunogenicity', pd.DataFrame())
        
        # 生成HTML内容
        html_content = self._generate_html(
            candidates, consensus, predictions, immunogenicity, species,
            epitope_candidates=epitope_candidates
        )
        
        # 写入文件
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"已生成HTML报告: {output_path}")
        
        return str(output_path)
    
    def _generate_html(self,
                       candidates: pd.DataFrame,
                       consensus: pd.DataFrame,
                       predictions: dict,
                       immunogenicity: pd.DataFrame,
                       species: Optional[str] = None,
                       epitope_candidates: pd.DataFrame = None) -> str:
        """生成HTML内容"""
        
        # 计算统计信息
        stats = self._calculate_stats(candidates, consensus, predictions)
        
        # 构建species显示信息
        species_display = ""
        if species:
            species_display = f'<p>目标细菌种类: <strong>{species}</strong></p>'
        
        # 生成图表
        if HAS_PLOTLY:
            charts = self._generate_plotly_charts(candidates, consensus)
        else:
            charts = self._generate_simple_charts(candidates, consensus)
        
        # 生成候选表格
        candidates_table = self._generate_candidates_table(candidates)
        
        # 生成表位级别候选表格
        epitope_table = ""
        if epitope_candidates is not None and not epitope_candidates.empty:
            epitope_table = self._generate_epitope_candidates_table(epitope_candidates)
        
        # 组装HTML
        html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BacterialAntigenFinder - 分析报告</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        header h1 {{
            font-size: 2.5rem;
            margin-bottom: 10px;
        }}
        header p {{
            opacity: 0.9;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .stat-card .value {{
            font-size: 2.5rem;
            font-weight: bold;
            color: #667eea;
        }}
        .stat-card .label {{
            color: #666;
            font-size: 0.9rem;
            margin-top: 5px;
        }}
        .section {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }}
        .section h2 {{
            color: #333;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #667eea;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #333;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .recommendation-high {{
            color: #28a745;
            font-weight: bold;
        }}
        .recommendation-medium {{
            color: #ffc107;
            font-weight: bold;
        }}
        .recommendation-low {{
            color: #dc3545;
        }}
        .chart-container {{
            margin: 20px 0;
        }}
        .footer {{
            text-align: center;
            color: #666;
            padding: 20px;
            font-size: 0.9rem;
        }}
        .epitope-sequence {{
            font-family: 'Courier New', monospace;
            background: #f0f0f0;
            padding: 2px 6px;
            border-radius: 3px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>BacterialAntigenFinder</h1>
            <p>细菌抗原AI智能筛选平台 - 分析报告</p>
            {species_display}
            <p>生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </header>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="value">{stats['total_proteins']}</div>
                <div class="label">分析蛋白质数</div>
            </div>
            <div class="stat-card">
                <div class="value">{stats['total_residues']}</div>
                <div class="label">总残基数</div>
            </div>
            <div class="stat-card">
                <div class="value">{stats['epitope_residues']}</div>
                <div class="label">表位残基数</div>
            </div>
            <div class="stat-card">
                <div class="value">{stats['high_priority']}</div>
                <div class="label">高优先级候选</div>
            </div>
            {'<div class="stat-card"><div class="value" style="font-size:1.2rem;">' + species + '</div><div class="label">细菌种类</div></div>' if species else ''}
        </div>
        
        <div class="section">
            <h2>Top 候选抗原（蛋白质级别）</h2>
            {candidates_table}
        </div>
        
        {'<div class="section"><h2>Top 候选抗原表位（表位级别）</h2>' + epitope_table + '</div>' if epitope_table else ''}
        
        <div class="section">
            <h2>预测分析</h2>
            <div class="chart-container">
                {charts}
            </div>
        </div>
        
        <div class="section">
            <h2>预测器统计</h2>
            {self._generate_predictor_stats(predictions)}
        </div>
        
        <footer class="footer">
            <p>Generated by BacterialAntigenFinder v1.0.0</p>
            <p>整合 BepiPred-3.0, DiscoTope-3.0, GraphBepi, EpiGraph, Vaxign-ML</p>
        </footer>
    </div>
</body>
</html>"""
        
        return html
    
    def _calculate_stats(self, candidates: pd.DataFrame, 
                          consensus: pd.DataFrame,
                          predictions: dict) -> dict:
        """计算统计信息"""
        stats = {
            'total_proteins': 0,
            'total_residues': 0,
            'epitope_residues': 0,
            'high_priority': 0,
            'medium_priority': 0,
            'low_priority': 0
        }
        
        if not consensus.empty:
            stats['total_proteins'] = consensus['protein_id'].nunique()
            stats['total_residues'] = len(consensus)
            stats['epitope_residues'] = consensus['is_consensus_epitope'].sum()
        
        if not candidates.empty and 'recommendation' in candidates.columns:
            stats['high_priority'] = (candidates['recommendation'] == 'HIGH').sum()
            stats['medium_priority'] = (candidates['recommendation'] == 'MEDIUM').sum()
            stats['low_priority'] = (candidates['recommendation'] == 'LOW').sum()
        
        return stats
    
    def _generate_candidates_table(self, candidates: pd.DataFrame) -> str:
        """生成候选抗原表格"""
        if candidates.empty:
            return "<p>暂无候选数据</p>"
        
        # 只显示前20个候选
        display_df = candidates.head(20)
        
        rows = []
        for _, row in display_df.iterrows():
            rec_class = f"recommendation-{row.get('recommendation', 'low').lower()}"
            
            rows.append(f"""
            <tr>
                <td>{row.get('rank', '-')}</td>
                <td><strong>{row.get('protein_id', '-')}</strong></td>
                <td>{row.get('residue_range', '-')}</td>
                <td><span class="epitope-sequence">{str(row.get('epitope_sequence', '-'))[:20]}...</span></td>
                <td>{row.get('avg_consensus_score', 0):.3f}</td>
                <td>{row.get('protegenicity_score', 0):.1f}</td>
                <td>{row.get('subcellular_location', '-')}</td>
                <td class="{rec_class}">{row.get('recommendation', '-')}</td>
            </tr>
            """)
        
        table = f"""
        <table>
            <thead>
                <tr>
                    <th>排名</th>
                    <th>蛋白质ID</th>
                    <th>最佳表位区域</th>
                    <th>表位序列</th>
                    <th>共识分数</th>
                    <th>保护性分数</th>
                    <th>亚细胞定位</th>
                    <th>推荐等级</th>
                </tr>
            </thead>
            <tbody>
                {''.join(rows)}
            </tbody>
        </table>
        """
        
        if len(candidates) > 20:
            table += f"<p style='margin-top:10px;color:#666;'>显示前20个候选，共{len(candidates)}个</p>"
        
        return table
    
    def _generate_epitope_candidates_table(self, epitope_candidates: pd.DataFrame) -> str:
        """生成表位级别候选抗原表格"""
        if epitope_candidates.empty:
            return "<p>暂无表位候选数据</p>"
        
        # 只显示前20个候选表位
        display_df = epitope_candidates.head(20)
        
        rows = []
        for _, row in display_df.iterrows():
            rec_class = f"recommendation-{row.get('recommendation', 'low').lower()}"
            
            rows.append(f"""
            <tr>
                <td>{row.get('rank', '-')}</td>
                <td><strong>{row.get('protein_id', '-')}</strong></td>
                <td>{row.get('epitope_start', '-')}-{row.get('epitope_end', '-')}</td>
                <td><span class="epitope-sequence">{str(row.get('epitope_sequence', '-'))[:25]}</span></td>
                <td>{row.get('epitope_length', '-')}</td>
                <td>{row.get('composite_score', 0):.4f}</td>
                <td>{row.get('avg_consensus_score', 0):.3f}</td>
                <td>{row.get('protegenicity_score', 0):.1f}</td>
                <td>{row.get('antigenicity_score', 0):.3f}</td>
                <td class="{rec_class}">{row.get('recommendation', '-')}</td>
            </tr>
            """)
        
        table = f"""
        <table>
            <thead>
                <tr>
                    <th>排名</th>
                    <th>蛋白质ID</th>
                    <th>表位位置</th>
                    <th>表位序列</th>
                    <th>长度</th>
                    <th>综合评分</th>
                    <th>共识分数</th>
                    <th>保护性分数</th>
                    <th>抗原性分数</th>
                    <th>推荐等级</th>
                </tr>
            </thead>
            <tbody>
                {''.join(rows)}
            </tbody>
        </table>
        """
        
        if len(epitope_candidates) > 20:
            table += f"<p style='margin-top:10px;color:#666;'>显示前20个候选表位，共{len(epitope_candidates)}个</p>"
        
        return table
    
    def _generate_predictor_stats(self, predictions: dict) -> str:
        """生成预测器统计表格"""
        if not predictions:
            return "<p>暂无预测器统计数据</p>"
        
        rows = []
        for name, result in predictions.items():
            if hasattr(result, 'predictions'):
                total = len(result.predictions)
                epitopes = sum(1 for p in result.predictions if p.is_epitope)
                ratio = epitopes / total * 100 if total > 0 else 0
                
                rows.append(f"""
                <tr>
                    <td><strong>{name}</strong></td>
                    <td>{total}</td>
                    <td>{epitopes}</td>
                    <td>{ratio:.1f}%</td>
                </tr>
                """)
        
        return f"""
        <table>
            <thead>
                <tr>
                    <th>预测器</th>
                    <th>总残基</th>
                    <th>表位残基</th>
                    <th>表位比例</th>
                </tr>
            </thead>
            <tbody>
                {''.join(rows)}
            </tbody>
        </table>
        """
    
    def _generate_plotly_charts(self, candidates: pd.DataFrame,
                                 consensus: pd.DataFrame) -> str:
        """生成Plotly图表"""
        if not HAS_PLOTLY or candidates.empty:
            return ""
        
        charts_html = []
        
        # 推荐分布饼图
        if 'recommendation' in candidates.columns:
            rec_counts = candidates['recommendation'].value_counts()
            fig = go.Figure(data=[go.Pie(
                labels=rec_counts.index.tolist(),
                values=rec_counts.values.tolist(),
                hole=0.4,
                marker_colors=['#28a745', '#ffc107', '#dc3545']
            )])
            fig.update_layout(
                title='候选抗原推荐分布',
                height=300,
                margin=dict(l=20, r=20, t=40, b=20)
            )
            charts_html.append(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        
        # 共识分数分布
        if not consensus.empty and 'consensus_score' in consensus.columns:
            fig = px.histogram(
                consensus, 
                x='consensus_score',
                nbins=50,
                title='共识分数分布'
            )
            fig.update_layout(
                height=300,
                margin=dict(l=20, r=20, t=40, b=20)
            )
            charts_html.append(fig.to_html(full_html=False, include_plotlyjs=False))
        
        return ''.join(charts_html)
    
    def _generate_simple_charts(self, candidates: pd.DataFrame,
                                 consensus: pd.DataFrame) -> str:
        """生成简单的HTML图表（无Plotly时使用）"""
        if candidates.empty:
            return ""
        
        # 使用CSS生成简单的条形图
        if 'recommendation' not in candidates.columns:
            return ""
        
        rec_counts = candidates['recommendation'].value_counts()
        total = rec_counts.sum()
        
        bars = []
        colors = {'HIGH': '#28a745', 'MEDIUM': '#ffc107', 'LOW': '#dc3545'}
        
        for rec, count in rec_counts.items():
            width = count / total * 100
            color = colors.get(rec, '#666')
            bars.append(f"""
            <div style="margin:5px 0;">
                <span style="display:inline-block;width:80px;">{rec}</span>
                <div style="display:inline-block;width:{width}%;background:{color};
                            height:20px;border-radius:3px;"></div>
                <span style="margin-left:10px;">{count}</span>
            </div>
            """)
        
        return f"""
        <div style="margin:20px 0;">
            <h3>候选抗原推荐分布</h3>
            {''.join(bars)}
        </div>
        """
