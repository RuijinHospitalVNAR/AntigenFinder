#!/usr/bin/env python3
"""
BacterialAntigenFinder 输出验证脚本

验证Pipeline运行后的输出文件是否完整和正确。

用法:
    python validate_output.py --output_dir results/
    python validate_output.py --output_dir results/ --species Pseudomonas_aeruginosa
    python validate_output.py --output_dir results/ --min_candidates 5
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional

import pandas as pd


# 颜色输出
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    NC = '\033[0m'


def print_pass(msg: str):
    print(f"  {Colors.GREEN}✓{Colors.NC} {msg}")


def print_fail(msg: str):
    print(f"  {Colors.RED}✗{Colors.NC} {msg}")


def print_warn(msg: str):
    print(f"  {Colors.YELLOW}⚠{Colors.NC} {msg}")


def print_info(msg: str):
    print(f"  {Colors.BLUE}ℹ{Colors.NC} {msg}")


def validate_antigen_candidates(output_dir: str, min_candidates: int = 5,
                                 expected_species: Optional[str] = None) -> dict:
    """
    验证antigen_candidates.csv文件

    Args:
        output_dir: 输出目录
        min_candidates: 最小候选数量
        expected_species: 预期的物种名称

    Returns:
        验证结果字典
    """
    result = {'passed': True, 'errors': [], 'warnings': [], 'info': []}
    csv_path = os.path.join(output_dir, 'antigen_candidates.csv')

    # 检查文件存在
    if not os.path.exists(csv_path):
        result['passed'] = False
        result['errors'].append('antigen_candidates.csv 不存在')
        return result

    print_pass('antigen_candidates.csv 存在')

    # 读取CSV
    try:
        # 跳过注释行
        df = pd.read_csv(csv_path, comment='#', encoding='utf-8-sig')
    except Exception as e:
        result['passed'] = False
        result['errors'].append(f'读取antigen_candidates.csv失败: {e}')
        return result

    # 检查候选数量
    num_candidates = len(df)
    if num_candidates < min_candidates:
        result['passed'] = False
        result['errors'].append(
            f'候选数量不足: {num_candidates} < {min_candidates}'
        )
        print_fail(f'候选数量不足: {num_candidates} < {min_candidates}')
    else:
        print_pass(f'候选数量: {num_candidates} (>= {min_candidates})')
    result['info'].append(f'num_candidates: {num_candidates}')

    # 检查必要列
    required_columns = [
        'Protein_ID', 'Recommendation', 'Composite_Score',
        'Avg_Consensus_Score', 'Rank'
    ]
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        # 尝试小写列名匹配
        alt_required = [col.lower() for col in required_columns]
        actual_cols_lower = {col.lower(): col for col in df.columns}
        still_missing = []
        for i, alt_col in enumerate(alt_required):
            if alt_col in actual_cols_lower:
                # 找到对应的大写列名
                pass
            else:
                still_missing.append(required_columns[i])
        if still_missing:
            result['passed'] = False
            result['errors'].append(f'缺少必要列: {still_missing}')
            print_fail(f'缺少必要列: {still_missing}')
        else:
            print_pass('所有必要列存在')
    else:
        print_pass('所有必要列存在')

    # 检查Composite_Score范围
    score_col = None
    for col_name in ['Composite_Score', 'composite_score']:
        if col_name in df.columns:
            score_col = col_name
            break

    if score_col:
        scores = df[score_col].dropna()
        if len(scores) > 0:
            min_score = scores.min()
            max_score = scores.max()
            if min_score < 0 or max_score > 1:
                result['passed'] = False
                result['errors'].append(
                    f'Composite_Score超出有效范围: [{min_score:.4f}, {max_score:.4f}]'
                )
                print_fail(f'Composite_Score超出有效范围: [{min_score:.4f}, {max_score:.4f}]')
            else:
                print_pass(f'Composite_Score范围有效: [{min_score:.4f}, {max_score:.4f}]')
            result['info'].append(f'composite_score_range: [{min_score:.4f}, {max_score:.4f}]')

    # 检查Recommendation值
    rec_col = None
    for col_name in ['Recommendation', 'recommendation']:
        if col_name in df.columns:
            rec_col = col_name
            break

    if rec_col:
        valid_recommendations = {'HIGH', 'MEDIUM', 'LOW'}
        actual_recs = set(df[rec_col].dropna().unique())
        invalid_recs = actual_recs - valid_recommendations
        if invalid_recs:
            result['passed'] = False
            result['errors'].append(f'无效的推荐等级: {invalid_recs}')
            print_fail(f'无效的推荐等级: {invalid_recs}')
        else:
            rec_counts = df[rec_col].value_counts().to_dict()
            print_pass(f'推荐等级有效: {rec_counts}')
            result['info'].append(f'recommendation_counts: {rec_counts}')

    # 检查Species列（如果预期）
    if expected_species:
        species_col = None
        for col_name in ['Species', 'species']:
            if col_name in df.columns:
                species_col = col_name
                break

        if species_col:
            species_values = df[species_col].dropna().unique()
            if expected_species in species_values:
                print_pass(f'Species列包含: {expected_species}')
            else:
                result['warnings'].append(
                    f'Species列不包含预期值 {expected_species}: {species_values}'
                )
                print_warn(f'Species列不包含预期值 {expected_species}: {species_values}')
        else:
            # 检查文件头注释
            try:
                with open(csv_path, 'r', encoding='utf-8-sig') as f:
                    first_line = f.readline().strip()
                if expected_species in first_line:
                    print_pass(f'文件头注释包含species: {expected_species}')
                else:
                    result['warnings'].append(
                        f'CSV文件中未找到species信息: {expected_species}'
                    )
                    print_warn(f'CSV文件中未找到species信息: {expected_species}')
            except Exception:
                pass

    return result


def validate_epitope_candidates(output_dir: str) -> dict:
    """
    验证epitope_candidates.csv文件

    Args:
        output_dir: 输出目录

    Returns:
        验证结果字典
    """
    result = {'passed': True, 'errors': [], 'warnings': [], 'info': []}

    # 检查主目录
    csv_path = os.path.join(output_dir, 'epitope_candidates.csv')
    if os.path.exists(csv_path):
        print_pass('epitope_candidates.csv 存在')
        try:
            df = pd.read_csv(csv_path, encoding='utf-8-sig')
            num_epitopes = len(df)
            print_pass(f'表位数量: {num_epitopes}')
            result['info'].append(f'num_epitopes: {num_epitopes}')

            # 检查必要列
            required_cols = ['protein_id', 'epitope_start', 'epitope_end',
                             'composite_score', 'recommendation']
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                result['warnings'].append(f'epitope_candidates.csv缺少列: {missing}')
                print_warn(f'缺少列: {missing}')
            else:
                print_pass('epitope_candidates.csv 所有必要列存在')

            # 检查推荐等级
            if 'recommendation' in df.columns:
                valid_recs = {'HIGH', 'MEDIUM', 'LOW'}
                actual_recs = set(df['recommendation'].dropna().unique())
                invalid = actual_recs - valid_recs
                if invalid:
                    result['passed'] = False
                    result['errors'].append(f'表位推荐等级无效: {invalid}')
                    print_fail(f'表位推荐等级无效: {invalid}')
                else:
                    print_pass('表位推荐等级有效')

        except Exception as e:
            result['warnings'].append(f'读取epitope_candidates.csv失败: {e}')
            print_warn(f'读取epitope_candidates.csv失败: {e}')
    else:
        # 检查detailed_results目录
        detailed_path = os.path.join(output_dir, 'detailed_results')
        if os.path.isdir(detailed_path):
            print_warn('epitope_candidates.csv 不存在于主目录，检查detailed_results/')
            result['warnings'].append('epitope_candidates.csv 不存在于主输出目录')
        else:
            result['warnings'].append('epitope_candidates.csv 不存在')
            print_warn('epitope_candidates.csv 不存在（可能是Pipeline未生成表位级别输出）')

    return result


def validate_run_metadata(output_dir: str,
                           expected_species: Optional[str] = None) -> dict:
    """
    验证run_metadata.json文件

    Args:
        output_dir: 输出目录
        expected_species: 预期的物种名称

    Returns:
        验证结果字典
    """
    result = {'passed': True, 'errors': [], 'warnings': [], 'info': []}
    json_path = os.path.join(output_dir, 'run_metadata.json')

    if not os.path.exists(json_path):
        result['passed'] = False
        result['errors'].append('run_metadata.json 不存在')
        return result

    print_pass('run_metadata.json 存在')

    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            metadata = json.load(f)
    except Exception as e:
        result['passed'] = False
        result['errors'].append(f'读取run_metadata.json失败: {e}')
        return result

    # 检查必要字段
    required_fields = ['version', 'species', 'input', 'output']
    for field in required_fields:
        if field not in metadata:
            result['warnings'].append(f'run_metadata.json缺少字段: {field}')
            print_warn(f'缺少字段: {field}')
        else:
            print_pass(f'包含字段: {field}')

    # 检查species信息
    species = metadata.get('species')
    if species:
        print_pass(f'species: {species}')
        result['info'].append(f'species: {species}')

        if expected_species and species != expected_species:
            result['passed'] = False
            result['errors'].append(
                f'species不匹配: 期望 "{expected_species}"，实际 "{species}"'
            )
            print_fail(f'species不匹配: 期望 "{expected_species}"，实际 "{species}"')
        elif expected_species:
            print_pass(f'species匹配: {expected_species}')
    else:
        result['warnings'].append('run_metadata.json中species为空')
        print_warn('species字段为空')

    # 检查版本
    version = metadata.get('version')
    if version:
        print_info(f'Pipeline版本: {version}')
        result['info'].append(f'version: {version}')

    # 检查运行时间
    elapsed = metadata.get('elapsed_seconds')
    if elapsed:
        print_info(f'运行耗时: {elapsed}秒')
        result['info'].append(f'elapsed_seconds: {elapsed}')

    # 检查输入信息
    input_info = metadata.get('input', {})
    if input_info:
        num_seq = input_info.get('num_sequences', 'N/A')
        print_info(f'输入序列数: {num_seq}')
        result['info'].append(f'num_sequences: {num_seq}')

    # 检查输出信息
    output_info = metadata.get('output', {})
    if output_info:
        num_cand = output_info.get('num_candidates', 'N/A')
        print_info(f'候选数量: {num_cand}')
        result['info'].append(f'num_candidates: {num_cand}')

    return result


def validate_html_report(output_dir: str) -> dict:
    """
    验证HTML报告文件

    Args:
        output_dir: 输出目录

    Returns:
        验证结果字典
    """
    result = {'passed': True, 'errors': [], 'warnings': [], 'info': []}
    html_path = os.path.join(output_dir, 'analysis_report.html')

    if not os.path.exists(html_path):
        result['warnings'].append('analysis_report.html 不存在')
        print_warn('analysis_report.html 不存在')
        return result

    print_pass('analysis_report.html 存在')

    try:
        with open(html_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # 检查基本HTML结构
        if '<html' in content and '</html>' in content:
            print_pass('HTML结构完整')
        else:
            result['warnings'].append('HTML结构可能不完整')
            print_warn('HTML结构可能不完整')

        # 检查关键内容
        if 'BacterialAntigenFinder' in content:
            print_pass('包含标题')
        if '候选抗原' in content or 'candidate' in content.lower():
            print_pass('包含候选数据')

        file_size = os.path.getsize(html_path)
        print_info(f'HTML文件大小: {file_size / 1024:.1f} KB')
        result['info'].append(f'html_size_kb: {file_size / 1024:.1f}')

    except Exception as e:
        result['warnings'].append(f'读取HTML报告失败: {e}')
        print_warn(f'读取HTML报告失败: {e}')

    return result


def validate_detailed_results(output_dir: str) -> dict:
    """
    验证详细结果目录

    Args:
        output_dir: 输出目录

    Returns:
        验证结果字典
    """
    result = {'passed': True, 'errors': [], 'warnings': [], 'info': []}
    detailed_dir = os.path.join(output_dir, 'detailed_results')

    if not os.path.isdir(detailed_dir):
        result['warnings'].append('detailed_results/ 目录不存在')
        print_warn('detailed_results/ 目录不存在')
        return result

    print_pass('detailed_results/ 目录存在')

    # 检查共识评分CSV
    consensus_path = os.path.join(detailed_dir, 'consensus_scores.csv')
    if os.path.exists(consensus_path):
        print_pass('consensus_scores.csv 存在')
    else:
        result['warnings'].append('consensus_scores.csv 不存在')
        print_warn('consensus_scores.csv 不存在')

    # 检查免疫原性CSV
    immuno_path = os.path.join(detailed_dir, 'immunogenicity_scores.csv')
    if os.path.exists(immuno_path):
        print_pass('immunogenicity_scores.csv 存在')
    else:
        result['warnings'].append('immunogenicity_scores.csv 不存在')
        print_warn('immunogenicity_scores.csv 不存在')

    # 列出目录内容
    files = os.listdir(detailed_dir)
    print_info(f'detailed_results/ 包含 {len(files)} 个文件/目录')
    result['info'].append(f'detailed_files: {len(files)}')

    return result


def run_validation(output_dir: str, min_candidates: int = 5,
                    expected_species: Optional[str] = None) -> bool:
    """
    运行完整的输出验证

    Args:
        output_dir: 输出目录
        min_candidates: 最小候选数量
        expected_species: 预期的物种名称

    Returns:
        验证是否通过
    """
    print(f"\n{'='*60}")
    print("BacterialAntigenFinder 输出验证")
    print(f"{'='*60}")
    print(f"输出目录: {output_dir}")
    if expected_species:
        print(f"预期物种: {expected_species}")
    print(f"最小候选数: {min_candidates}")
    print()

    # 检查输出目录存在
    if not os.path.isdir(output_dir):
        print_fail(f'输出目录不存在: {output_dir}')
        return False

    all_passed = True
    all_info = []
    all_errors = []
    all_warnings = []

    # 1. 验证antigen_candidates.csv
    print("--- 验证 antigen_candidates.csv ---")
    r = validate_antigen_candidates(output_dir, min_candidates, expected_species)
    all_passed = all_passed and r['passed']
    all_errors.extend(r['errors'])
    all_warnings.extend(r['warnings'])
    all_info.extend(r['info'])
    print()

    # 2. 验证epitope_candidates.csv
    print("--- 验证 epitope_candidates.csv ---")
    r = validate_epitope_candidates(output_dir)
    all_passed = all_passed and r['passed']
    all_errors.extend(r['errors'])
    all_warnings.extend(r['warnings'])
    all_info.extend(r['info'])
    print()

    # 3. 验证run_metadata.json
    print("--- 验证 run_metadata.json ---")
    r = validate_run_metadata(output_dir, expected_species)
    all_passed = all_passed and r['passed']
    all_errors.extend(r['errors'])
    all_warnings.extend(r['warnings'])
    all_info.extend(r['info'])
    print()

    # 4. 验证HTML报告
    print("--- 验证 analysis_report.html ---")
    r = validate_html_report(output_dir)
    # HTML报告缺失不是致命错误
    all_errors.extend(r['errors'])
    all_warnings.extend(r['warnings'])
    all_info.extend(r['info'])
    print()

    # 5. 验证详细结果
    print("--- 验证 detailed_results/ ---")
    r = validate_detailed_results(output_dir)
    all_errors.extend(r['errors'])
    all_warnings.extend(r['warnings'])
    all_info.extend(r['info'])
    print()

    # 打印汇总
    print(f"{'='*60}")
    print("验证汇总")
    print(f"{'='*60}")

    if all_passed:
        print(f"\n  {Colors.GREEN}验证通过!{Colors.NC}")
    else:
        print(f"\n  {Colors.RED}验证失败!{Colors.NC}")

    if all_errors:
        print(f"\n  错误 ({len(all_errors)}):")
        for err in all_errors:
            print(f"    {Colors.RED}- {err}{Colors.NC}")

    if all_warnings:
        print(f"\n  警告 ({len(all_warnings)}):")
        for warn in all_warnings:
            print(f"    {Colors.YELLOW}- {warn}{Colors.NC}")

    if all_info:
        print(f"\n  信息:")
        for info in all_info:
            print(f"    {Colors.BLUE}- {info}{Colors.NC}")

    print()
    return all_passed


def main():
    parser = argparse.ArgumentParser(
        description='BacterialAntigenFinder 输出验证工具'
    )
    parser.add_argument(
        '--output_dir', '-o',
        type=str,
        required=True,
        help='Pipeline输出目录'
    )
    parser.add_argument(
        '--species', '-s',
        type=str,
        default=None,
        help='预期的细菌种类名称 (如 Pseudomonas_aeruginosa)'
    )
    parser.add_argument(
        '--min_candidates',
        type=int,
        default=5,
        help='最小候选数量 (默认: 5)'
    )

    args = parser.parse_args()

    passed = run_validation(
        output_dir=args.output_dir,
        min_candidates=args.min_candidates,
        expected_species=args.species
    )

    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
