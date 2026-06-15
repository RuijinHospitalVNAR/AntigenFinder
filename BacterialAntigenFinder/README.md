# BacterialAntigenFinder

细菌抗原AI智能筛选平台 - 整合多种B细胞表位预测模型，通过共识评分和免疫原性评估，输出针对耐药细菌的抗原候选清单。

## 功能特点

- **多模型集成**: 整合 BepiPred-3.0、DiscoTope-3.0、GraphBepi、EpiGraph、Vaxign-ML 五种预测模型
- **共识评分**: 加权平均 + 多数投票机制，提高预测准确性
- **免疫原性评估**: 基于 Vaxign-ML 的保护性抗原评估
- **灵活配置**: 支持自定义模型权重、阈值和输出格式
- **Docker支持**: 提供完整的 Docker 部署方案

## 安装方式

### 方式一: Conda环境安装

```bash
# 1. 克隆项目
git clone [repo_url]
cd BacterialAntigenFinder

# 2. 运行环境安装脚本
chmod +x scripts/setup_envs.sh
./scripts/setup_envs.sh

# 3. 激活主控环境
conda activate master_env
```

### 方式二: Docker安装

```bash
# 1. 构建镜像
docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .

# 2. 运行容器
docker run --gpus all \
    -v $(pwd)/data:/data \
    -v $(pwd)/results:/results \
    bacterial-antigen-finder:latest \
    --fasta /data/antigens.fasta \
    --pdb_dir /data/structures/ \
    --organism_type gram- \
    --output_dir /results/
```

## 使用方法

### 基本用法

```bash
python main.py \
    --fasta input_sequences.fasta \
    --pdb_dir structures/ \
    --organism_type gram- \
    --output_dir results/
```

### 完整参数

```bash
python main.py \
    --fasta input_sequences.fasta \
    --pdb_dir structures/ \
    --organism_type gram- \
    --output_dir results/ \
    --models bepipred,discotope,graphbepi,epigraph,vaxignml \
    --consensus_threshold 0.5 \
    --top_candidates 50 \
    --output_format both
```

### 参数说明

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--fasta, -f` | str | 必需 | FASTA格式的蛋白质序列文件 |
| `--pdb_dir, -p` | str | 必需 | PDB结构文件目录 |
| `--organism_type, -t` | str | gram- | 生物类型 (gram+, gram-, virus) |
| `--models, -m` | str | 全部 | 启用的模型，逗号分隔 |
| `--output_dir, -o` | str | ./results | 输出目录 |
| `--consensus_threshold` | float | 0.5 | 共识评分阈值 |
| `--min_votes` | int | 2 | 最小投票数 |
| `--top_candidates` | int | 50 | 输出前N个候选 |
| `--output_format` | str | both | 输出格式 (csv, html, both) |
| `--cpu_only` | flag | False | 仅使用CPU |
| `--config, -c` | str | None | 自定义配置文件路径 |

## 输入格式

### FASTA序列文件

```fasta
>OmpA_ECOLI Outer membrane protein A
MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNG...
>FimH_ECOLI Type 1 fimbrial adhesin
MKRVITLFAVLLMGWSVNAWSFACKTANGTAIPIGGGSANVYVNLAPVVNVG...
```

### PDB结构文件

- 支持标准PDB格式
- 文件名应与FASTA中的序列ID对应（如 `OmpA_ECOLI.pdb`）
- 支持实验解析结构和AlphaFold预测结构

## 输出格式

### 抗原候选清单 (CSV)

```csv
protein_id,residue_range,epitope_sequence,bepipred_score,discotope_score,graphbepi_score,epigraph_score,consensus_score,vote_count,vaxignml_protegenicity,immunogenicity_rank,recommendation
OmpA_ECOLI,45-60,KGFNKTGDTGVLQ,0.82,0.75,0.88,0.71,0.79,4,85.2,1,HIGH
FimH_ECOLI,120-135,NNPVTGQGT...,0.65,0.80,0.72,0.68,0.71,4,78.5,2,HIGH
```

### HTML分析报告

- 交互式可视化图表
- 各模型预测结果对比
- 候选抗原详细信息

## 模型说明

| 模型 | 类型 | 输入 | 说明 |
|------|------|------|------|
| BepiPred-3.0 | 序列基 | FASTA | 基于ESM-2的B细胞表位预测 |
| DiscoTope-3.0 | 结构基 | PDB | 基于ESM-IF1的构象表位预测 |
| GraphBepi | 混合 | PDB | 图神经网络表位预测 |
| EpiGraph | 结构基 | PDB | 图注意力网络表位预测 |
| Vaxign-ML | 序列基 | FASTA | 保护性抗原预测和免疫原性评估 |

## 项目结构

```
BacterialAntigenFinder/
├── main.py                      # 主入口CLI
├── config/
│   └── default_config.yaml      # 默认配置
├── src/
│   ├── preprocessor/            # 预处理模块
│   ├── predictors/              # 预测器包装器
│   ├── aggregator/              # 结果聚合模块
│   └── reporter/                # 报告生成模块
├── envs/                        # Conda环境配置
├── docker/                      # Docker部署文件
├── scripts/                     # 安装和运行脚本
├── tests/                       # 测试用例
└── example_data/                # 示例数据
```

## 引用

如果您使用了本工具，请引用相关模型的论文：

- BepiPred-3.0: Clifford et al., Protein Science (2022)
- DiscoTope-3.0: Høie et al., Frontiers in Immunology (2024)
- GraphBepi: Zeng et al., bioRxiv (2022)
- EpiGraph: Choi & Kim (2023)
- Vaxign-ML: Ong et al. (2020)

## 许可证

本项目采用 MIT 许可证。各子模型请遵循其原始许可证要求。

## 联系方式

如有问题或建议，请提交 Issue 或联系开发团队。
