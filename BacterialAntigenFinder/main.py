#!/usr/bin/env python3
"""
BacterialAntigenFinder - 细菌抗原AI智能筛选平台

主入口CLI程序，整合多种B细胞表位预测模型进行抗原筛选。

版本: 1.0.0
作者: BacterialAntigenFinder Team
"""

import argparse
import logging
import sys
import os
import time
import json
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml

# 添加src到路径
sys.path.insert(0, str(Path(__file__).parent))

from src.preprocessor import FastaParser, PdbValidator, DataMapper
from src.predictors import (
    BepipredWrapper,
    DiscotopeWrapper,
    GraphbepiWrapper,
    EpigraphWrapper,
    VaxignmlWrapper
)
from src.aggregator import ConsensusScorer, ImmunogenicityEvaluator, CandidateRanker
from src.reporter import CsvExporter, HtmlReporter

__version__ = "1.0.0"


def setup_logging(log_level: str, log_file: Optional[str] = None):
    """配置日志"""
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, encoding='utf-8'))
    
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def load_config(config_path: str) -> dict:
    """加载配置文件"""
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='BacterialAntigenFinder - 细菌抗原AI智能筛选平台',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 基本用法
  python main.py --fasta antigens.fasta --pdb_dir structures/ --output_dir results/

  # 指定生物类型和模型
  python main.py --fasta antigens.fasta --pdb_dir structures/ \\
      --organism_type gram- --models bepipred,discotope,graphbepi

  # 使用自定义配置
  python main.py --config my_config.yaml --fasta antigens.fasta
        """
    )
    
    # 输入参数
    input_group = parser.add_argument_group('输入参数')
    input_group.add_argument(
        '--fasta', '-f',
        type=str,
        required=True,
        help='FASTA格式的蛋白质序列文件'
    )
    input_group.add_argument(
        '--pdb_dir', '-p',
        type=str,
        required=True,
        help='PDB结构文件目录'
    )
    input_group.add_argument(
        '--organism_type', '-t',
        type=str,
        choices=['gram+', 'gram-', 'g+', 'g-', 'virus', 'v'],
        default='gram-',
        help='生物类型 (默认: gram-)'
    )
    input_group.add_argument(
        '--species', '-s',
        type=str,
        default=None,
        help='目标细菌种类 (如 Pseudomonas_aeruginosa, Klebsiella_pneumoniae)'
    )
    
    # 模型参数
    model_group = parser.add_argument_group('模型参数')
    model_group.add_argument(
        '--models', '-m',
        type=str,
        default='bepipred,discotope,graphbepi,epigraph,vaxignml',
        help='启用的模型，逗号分隔 (默认: 全部)'
    )
    model_group.add_argument(
        '--consensus_threshold',
        type=float,
        default=0.5,
        help='共识评分阈值 (默认: 0.5)'
    )
    model_group.add_argument(
        '--min_votes',
        type=int,
        default=2,
        help='最小投票数 (默认: 2)'
    )
    
    # 输出参数
    output_group = parser.add_argument_group('输出参数')
    output_group.add_argument(
        '--output_dir', '-o',
        type=str,
        default='./results',
        help='输出目录 (默认: ./results)'
    )
    output_group.add_argument(
        '--top_candidates',
        type=int,
        default=50,
        help='输出前N个候选 (默认: 50)'
    )
    output_group.add_argument(
        '--output_format',
        type=str,
        choices=['csv', 'html', 'both'],
        default='both',
        help='输出格式 (默认: both)'
    )
    
    # 运行参数
    runtime_group = parser.add_argument_group('运行参数')
    runtime_group.add_argument(
        '--config', '-c',
        type=str,
        default=None,
        help='配置文件路径'
    )
    runtime_group.add_argument(
        '--cpu_only',
        action='store_true',
        help='仅使用CPU'
    )
    runtime_group.add_argument(
        '--workers',
        type=int,
        default=4,
        help='并行工作线程数 (默认: 4)'
    )
    runtime_group.add_argument(
        '--continue_on_error',
        action='store_true',
        default=True,
        help='单个模型失败时继续运行'
    )
    runtime_group.add_argument(
        '--log_level',
        type=str,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='日志级别 (默认: INFO)'
    )
    runtime_group.add_argument(
        '--log_file',
        type=str,
        default=None,
        help='日志文件路径'
    )
    
    return parser.parse_args()


class AntigenFinderPipeline:
    """抗原筛选Pipeline主类"""
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        self.species = config.get('species', {}).get('name', None)
        
        # 初始化预处理器
        self.fasta_parser = FastaParser()
        self.pdb_validator = PdbValidator()
        self.data_mapper = DataMapper()
        
        # 初始化预测器
        self.predictors = {}
        self._init_predictors()
        
        # 初始化聚合器
        self.consensus_scorer = ConsensusScorer(
            weights=config.get('models', {}).get('weights', {}),
            threshold=config.get('consensus', {}).get('threshold', 0.5)
        )
        
        # 确定生物类型用于VaxiJen模型
        organism_type = config.get('organism_type', 'bacteria')
        # 规范化organism_type
        if organism_type in ['gram+', 'gram-', 'g+', 'g-']:
            organism_type = 'bacteria'
        elif organism_type in ['virus', 'v']:
            organism_type = 'virus'
        
        self.immunogenicity_evaluator = ImmunogenicityEvaluator(
            protegenicity_threshold=config.get('immunogenicity', {}).get('protegenicity_threshold', 50.0),
            antigenicity_threshold=config.get('immunogenicity', {}).get('antigenicity_threshold', 0.4),
            organism_type=organism_type,
            use_vaxijen=config.get('immunogenicity', {}).get('use_vaxijen', True)
        )
        
        self.candidate_ranker = CandidateRanker(
            top_n=config.get('output', {}).get('top_candidates', 50)
        )
        
        # 初始化报告器
        self.csv_exporter = CsvExporter()
        self.html_reporter = HtmlReporter()
    
    def _init_predictors(self):
        """初始化预测器"""
        enabled_models = self.config.get('models', {}).get('enabled', [])
        model_paths = self.config.get('model_paths', {})
        conda_envs = self.config.get('conda_envs', {})
        thresholds = self.config.get('models', {}).get('thresholds', {})
        use_gpu = self.config.get('runtime', {}).get('use_gpu', True)
        
        predictor_classes = {
            'bepipred': BepipredWrapper,
            'discotope': DiscotopeWrapper,
            'graphbepi': GraphbepiWrapper,
            'epigraph': EpigraphWrapper,
            'vaxignml': VaxignmlWrapper
        }
        
        for model_name in enabled_models:
            if model_name in predictor_classes:
                try:
                    self.predictors[model_name] = predictor_classes[model_name](
                        model_path=model_paths.get(model_name, ''),
                        env_name=conda_envs.get(model_name, f'{model_name}_env'),
                        threshold=thresholds.get(model_name, 0.5),
                        use_gpu=use_gpu
                    )
                    self.logger.info(f"已初始化预测器: {model_name}")
                except Exception as e:
                    self.logger.warning(f"初始化预测器 {model_name} 失败: {e}")
    
    def run(self, fasta_file: str, pdb_dir: str, organism_type: str, 
            output_dir: str, output_format: str = 'both') -> dict:
        """
        运行完整的抗原筛选流程
        
        Args:
            fasta_file: FASTA序列文件路径
            pdb_dir: PDB结构文件目录
            organism_type: 生物类型
            output_dir: 输出目录
            output_format: 输出格式
        
        Returns:
            包含候选抗原列表的字典
        """
        pipeline_start = time.time()
        
        self.logger.info("=" * 60)
        self.logger.info("BacterialAntigenFinder - 开始抗原筛选流程")
        self.logger.info(f"版本: {__version__}")
        self.logger.info(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("=" * 60)
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # Step 1: 预处理
        self.logger.info("Step 1: 数据预处理...")
        sequences = self.fasta_parser.parse(fasta_file)
        self.logger.info(f"  解析到 {len(sequences)} 个序列")
        
        structures = self.pdb_validator.validate_directory(pdb_dir)
        self.logger.info(f"  验证通过 {len(structures)} 个PDB结构")
        
        mapped_data = self.data_mapper.map(sequences, structures)
        self.logger.info(f"  成功映射 {len(mapped_data)} 个序列-结构对")
        
        # Step 2: 运行各预测器
        self.logger.info("Step 2: 运行表位预测模型...")
        predictions = {}
        continue_on_error = self.config.get('runtime', {}).get('continue_on_error', True)
        
        for name, predictor in self.predictors.items():
            try:
                start_time = time.time()
                self.logger.info(f"  运行 {name}...")
                
                # 根据预测器类型决定是否传入结构
                if predictor.input_type == 'sequence':
                    result = predictor.predict(
                        sequences=sequences,
                        structures=None,
                        organism_type=organism_type
                    )
                else:
                    result = predictor.predict(
                        sequences=sequences,
                        structures=structures,
                        organism_type=organism_type
                    )
                
                predictions[name] = result
                elapsed = time.time() - start_time
                pred_count = len(result.predictions) if hasattr(result, 'predictions') else 0
                epitope_count = sum(1 for p in result.predictions if p.is_epitope) if hasattr(result, 'predictions') else 0
                self.logger.info(f"  {name} 完成 ({elapsed:.1f}s)，{pred_count} 残基，{epitope_count} 表位")
            except Exception as e:
                self.logger.error(f"  {name} 预测失败: {e}")
                if not continue_on_error:
                    raise
        
        # Step 3: 共识评分
        self.logger.info("Step 3: 计算共识评分...")
        consensus_results = self.consensus_scorer.compute_consensus(predictions)
        self.logger.info(f"  共识评分完成，{len(consensus_results)} 个残基")
        
        # Step 4: 免疫原性评估
        self.logger.info("Step 4: 免疫原性评估...")
        if 'vaxignml' in predictions:
            immunogenicity_scores = self.immunogenicity_evaluator.evaluate(
                predictions['vaxignml'],
                sequences
            )
            self.logger.info(f"  免疫原性评估完成")
        else:
            immunogenicity_scores = None
            self.logger.warning("  跳过免疫原性评估（Vaxign-ML未运行）")
        
        # Step 5: 候选排序
        self.logger.info("Step 5: 候选抗原排序...")
        ranked_candidates = self.candidate_ranker.rank(
            consensus_results,
            immunogenicity_scores,
            sequences
        )
        self.logger.info(f"  蛋白质级别排序完成，Top {len(ranked_candidates)} 候选")
        
        # Step 5b: 表位级别排序
        self.logger.info("Step 5b: 表位级别排序...")
        ranked_epitopes = self.candidate_ranker.rank_epitopes(
            consensus_results,
            immunogenicity_scores,
            sequences
        )
        self.logger.info(f"  表位级别排序完成，Top {len(ranked_epitopes)} 候选表位")
        
        # Step 6: 输出结果
        self.logger.info("Step 6: 生成结果报告...")
        results = {
            'candidates': ranked_candidates,
            'epitope_candidates': ranked_epitopes,
            'consensus': consensus_results,
            'predictions': predictions,
            'immunogenicity': immunogenicity_scores
        }
        
        if output_format in ['csv', 'both']:
            csv_path = self.csv_exporter.export(
                ranked_candidates,
                os.path.join(output_dir, 'antigen_candidates.csv'),
                species=self.species
            )
            self.logger.info(f"  CSV报告: {csv_path}")
            
            # 导出表位级别候选清单
            if not ranked_epitopes.empty:
                epitope_csv_path = self.csv_exporter.export(
                    ranked_epitopes,
                    os.path.join(output_dir, 'epitope_candidates.csv'),
                    species=self.species
                )
                self.logger.info(f"  表位级别CSV: {epitope_csv_path}")
        
        if output_format in ['html', 'both']:
            html_path = self.html_reporter.generate(
                results,
                os.path.join(output_dir, 'analysis_report.html'),
                species=self.species
            )
            self.logger.info(f"  HTML报告: {html_path}")
        
        # 保存详细结果
        detailed_dir = os.path.join(output_dir, 'detailed_results')
        os.makedirs(detailed_dir, exist_ok=True)
        
        # 保存共识结果
        if not consensus_results.empty:
            consensus_results.to_csv(
                os.path.join(detailed_dir, 'consensus_scores.csv'),
                index=False
            )
        
        # 保存免疫原性结果
        if immunogenicity_scores is not None and not immunogenicity_scores.empty:
            immunogenicity_scores.to_csv(
                os.path.join(detailed_dir, 'immunogenicity_scores.csv'),
                index=False
            )
        
        # 保存运行元数据
        pipeline_elapsed = time.time() - pipeline_start
        metadata = {
            'version': __version__,
            'run_time': datetime.now().isoformat(),
            'elapsed_seconds': round(pipeline_elapsed, 2),
            'species': self.species,
            'input': {
                'fasta_file': fasta_file,
                'pdb_dir': pdb_dir,
                'organism_type': organism_type,
                'num_sequences': len(sequences),
                'num_structures': len(structures)
            },
            'models': list(self.predictors.keys()),
            'output': {
                'num_candidates': len(ranked_candidates),
                'num_epitope_candidates': len(ranked_epitopes),
                'output_dir': output_dir
            }
        }
        
        with open(os.path.join(output_dir, 'run_metadata.json'), 'w') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        
        self.logger.info("=" * 60)
        self.logger.info("抗原筛选流程完成!")
        self.logger.info(f"总耗时: {pipeline_elapsed:.1f} 秒")
        self.logger.info(f"结果保存至: {output_dir}")
        self.logger.info("=" * 60)
        
        return results


def main():
    """主函数"""
    args = parse_args()
    
    # 设置日志
    setup_logging(args.log_level, args.log_file)
    logger = logging.getLogger('main')
    
    # 加载配置
    if args.config:
        config = load_config(args.config)
    else:
        # 使用默认配置
        default_config_path = Path(__file__).parent / 'config' / 'default_config.yaml'
        if default_config_path.exists():
            config = load_config(str(default_config_path))
        else:
            config = {}
    
    # 命令行参数覆盖配置
    enabled_models = args.models.split(',')
    config.setdefault('models', {})['enabled'] = enabled_models
    config.setdefault('consensus', {})['threshold'] = args.consensus_threshold
    config.setdefault('consensus', {})['min_votes'] = args.min_votes
    config.setdefault('output', {})['top_candidates'] = args.top_candidates
    config.setdefault('runtime', {})['use_gpu'] = not args.cpu_only
    config.setdefault('runtime', {})['max_workers'] = args.workers
    config.setdefault('runtime', {})['continue_on_error'] = args.continue_on_error
    
    # 设置细菌种类
    if args.species:
        config.setdefault('species', {})['name'] = args.species
    
    # 规范化生物类型
    organism_type_map = {'g+': 'gram+', 'g-': 'gram-', 'v': 'virus'}
    organism_type = organism_type_map.get(args.organism_type, args.organism_type)
    
    try:
        # 创建并运行Pipeline
        pipeline = AntigenFinderPipeline(config)
        results = pipeline.run(
            fasta_file=args.fasta,
            pdb_dir=args.pdb_dir,
            organism_type=organism_type,
            output_dir=args.output_dir,
            output_format=args.output_format
        )
        
        # 打印摘要
        print(f"\n{'='*60}")
        print("抗原筛选完成!")
        print(f"{'='*60}")
        
        # 获取species显示信息
        species_name = args.species
        if species_name:
            predefined = config.get('species', {}).get('predefined', {})
            species_info = predefined.get(species_name, {})
            display_name = species_info.get('display_name', species_name)
            drug_resistance = species_info.get('drug_resistance', '')
            if drug_resistance:
                print(f"细菌种类: {display_name} ({species_name}) [{drug_resistance}]")
            else:
                print(f"细菌种类: {display_name} ({species_name})")
        
        # 统计序列数
        consensus_df = results.get('consensus')
        if consensus_df is not None and not consensus_df.empty:
            num_proteins = consensus_df['protein_id'].nunique()
            num_epitopes = consensus_df['is_consensus_epitope'].sum()
        else:
            num_proteins = 0
            num_epitopes = 0
        
        candidates_df = results.get('candidates')
        epitope_df = results.get('epitope_candidates')
        num_candidates = len(candidates_df) if candidates_df is not None else 0
        num_epitope_candidates = len(epitope_df) if epitope_df is not None else 0
        
        print(f"分析蛋白质数: {num_proteins}")
        print(f"共识表位残基: {num_epitopes}")
        print(f"候选抗原数: {num_candidates}")
        print(f"候选表位数: {num_epitope_candidates}")
        print(f"结果目录: {args.output_dir}")
        print(f"{'='*60}\n")
        
    except Exception as e:
        logger.error(f"Pipeline运行失败: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
