# AntigenFinder

**细菌抗原AI智能筛选平台** — 针对耐药细菌，整合多种B细胞表位预测模型，通过共识评分与免疫原性评估，智能输出抗原候选清单。

---

## 项目目标

耐药细菌（如铜绿假单胞菌、肺炎克雷伯菌等）的疫苗研发亟需高效的抗原筛选工具。本项目旨在构建一个**具备免疫原性评估功能的AI抗原智能筛选平台**，实现：

- 针对多种耐药细菌，自动化筛选高免疫原性抗原表位
- 整合多种AI预测模型，通过共识机制提升预测可靠性
- 输出带综合AI评分的候选抗原/表位清单，为疫苗设计提供决策支持

## 主要内容

### 支持的耐药细菌

| 细菌种类 | 中文名 | 耐药类型 | 革兰氏 |
|----------|--------|----------|--------|
| *Pseudomonas aeruginosa* | 铜绿假单胞菌 | MDR/XDR | G- |
| *Klebsiella pneumoniae* | 肺炎克雷伯菌 | ESBL/CRE | G- |
| *Acinetobacter baumannii* | 鲍曼不动杆菌 | XDR | G- |
| *Staphylococcus aureus* | 金黄色葡萄球菌 | MRSA | G+ |
| *Escherichia coli* | 大肠杆菌 | ESBL | G- |

### 核心功能

1. **多模型集成预测** — 整合4种B细胞表位预测模型
2. **共识评分** — 加权平均 + 多数投票；支持百分位排名归一化，降低单一模型偏差并解决不同模型分数尺度不一致问题
3. **免疫原性评估** — VaxiJen 2.0抗原性计算
4. **双粒度排序** — 蛋白质级别 & 表位级别候选排序
5. **AI综合评分** — 多维度加权评分，输出推荐等级（HIGH/MEDIUM/LOW）
6. **Docker部署** — 完整容器化方案，支持第三方复现

## 实现方式

### 集成模型

| 模型 | 类型 | 输入 | 核心算法 | 权重 |
|------|------|------|----------|------|
| BepiPred-3.0 | 序列基 | FASTA | ESM-2蛋白质语言模型 | 0.30 |
| DiscoTope-3.0 | 结构基 | PDB | ESM-IF1 + 结构特征 | 0.30 |
| GraphBepi | 结构基 | PDB | 图神经网络(GNN) | 0.20 |
| EpiGraph | 结构基 | PDB | 图注意力网络(GAT) | 0.20 |

### 评分算法

**共识评分**：对每个残基，按模型权重计算加权平均分，结合多数投票（≥2个模型同意）确定表位。由于不同模型输出的分数尺度差异较大（如 BepiPred 输出 0-1 概率，GraphBepi/EpiGraph 输出 GNN/GAT logit），默认启用 **百分位排名归一化**（`--score_normalize rank`），将每个模型内部分数映射到 `[0,1]` 后再聚合。

**免疫原性评估**：
- 抗原性分数（VaxiJen 2.0，基于氨基酸 z-scale 自相关系数）
- 亚细胞定位权重（外膜/分泌蛋白优先）
- 跨膜区域降权

**综合评分公式**：

蛋白质级别：`0.40×共识分 + 0.25×抗原性分 + 0.20×残基数量 + 0.15×区域数量`

表位级别：`0.35×共识分 + 0.25×抗原性分 + 0.20×表位长度 + 0.20×模型一致性`

### 主要流程

```
┌─────────────────────────────────────────────────────────────────┐
│                     AntigenFinder Pipeline                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  输入                                                            │
│  ├── FASTA 蛋白质序列文件                                        │
│  └── PDB 蛋白质结构文件                                          │
│       │                                                          │
│       ▼                                                          │
│  Step 1: 数据预处理                                              │
│  ├── FASTA序列解析 & 序列验证                                    │
│  ├── PDB结构验证                                                 │
│  └── 序列-结构映射                                               │
│       │                                                          │
│       ▼                                                          │
│  Step 2: 多模型并行预测                                          │
│  ├── BepiPred-3.0  ──→  序列基B细胞表位预测                     │
│  ├── DiscoTope-3.0  ──→  结构基构象表位预测                     │
│  ├── GraphBepi      ──→  GNN表位预测                            │
│  └── EpiGraph       ──→  GAT表位预测                            │
│       │                                                          │
│       ▼                                                          │
│  Step 3: 共识评分                                                │
│  ├── 加权平均 (权重: 0.30/0.30/0.20/0.20)                       │
│  ├── 多数投票 (≥2 models agree)                                 │
│  └── 共识表位识别                                                │
│       │                                                          │
│       ▼                                                          │
│  Step 4: 免疫原性评估                                            │
│  ├── VaxiJen 2.0 抗原性计算 (z-scale ACC)                       │
│  ├── 亚细胞定位分析                                              │
│  └── 推荐等级判定 (HIGH/MEDIUM/LOW)                              │
│       │                                                          │
│       ▼                                                          │
│  Step 5: 候选排序                                                │
│  ├── 蛋白质级别排序 → antigen_candidates.csv                     │
│  └── 表位级别排序   → epitope_candidates.csv                     │
│       │                                                          │
│       ▼                                                          │
│  Step 6: 报告输出                                                │
│  ├── CSV 候选清单 (蛋白质级 + 表位级)                            │
│  ├── HTML 分析报告 (可视化图表)                                  │
│  ├── 详细残基级预测结果                                          │
│  └── 运行元数据 (JSON)                                           │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## 安装

### 方式一：Conda环境安装（推荐）

```bash
# 1. 克隆项目
git clone https://github.com/RuijinHospitalVNAR/AntigenFinder.git
cd AntigenFinder/BacterialAntigenFinder-2

# 2. 运行环境安装脚本
chmod +x scripts/setup_envs.sh
./scripts/setup_envs.sh

# 3. 激活主控环境
conda activate master_env
```

### 方式二：Docker安装（推荐用于Linux服务器）

#### 前置条件

- Linux 系统（Ubuntu 20.04+ 推荐）
- Docker 20.10+ 和 Docker Compose v2
- NVIDIA GPU + 驱动 + [nvidia-container-toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)（GPU运行需要）
- 磁盘空间 ≥ 50GB（镜像约30GB + 模型文件）

#### Step 1: 克隆项目

```bash
git clone https://github.com/RuijinHospitalVNAR/AntigenFinder.git
cd AntigenFinder/BacterialAntigenFinder-2
```

#### Step 2: 准备预测模型

4个预测模型需要从 GitHub 克隆，放置在项目同级目录的 `models/` 下：

```bash
cd ..
mkdir -p models && cd models

# 克隆4个模型（BepiPred需可写，不加 --depth 1 以便缓存ESM编码）
git clone https://github.com/mabu-car/BepiPred-3.0.git BepiPred-3.0-main
git clone https://github.com/mabu-car/DiscoTope-3.0.git DiscoTope-3.0-master
git clone https://github.com/GraphBepi/GraphBepi.git GraphBepi-main
git clone https://github.com/GraphBepi/EpiGraph.git EpiGraph-main

cd ../BacterialAntigenFinder-2
```

最终目录结构：
```
AntigenFinder/
├── BacterialAntigenFinder-2/     # 项目代码
└── models/                       # 4个预测模型
    ├── BepiPred-3.0-main/
    ├── DiscoTope-3.0-master/
    ├── GraphBepi-main/
    └── EpiGraph-main/
```

#### Step 3: 准备PDB结构文件

DiscoTope/GraphBepi/EpiGraph 需要 PDB 结构文件。从 [AlphaFold DB](https://alphafold.ebi.ac.uk/) 下载对应蛋白结构：

```bash
mkdir -p example_data/structures

# 示例：下载铜绿假单胞菌 OprF (UniProt: P02443) 的AlphaFold预测结构
wget https://alphafold.ebi.ac.uk/files/AF-P02443-F1-model_v4.pdb \
    -O example_data/structures/OprF_PSEAE.pdb

# 为FASTA中每个蛋白下载对应PDB（文件名需与FASTA header中的protein_id一致）
```

> **提示**：如果暂时没有PDB文件，可以只运行序列基的 BepiPred 模型（见下方"仅序列模式"）。

#### Step 4: 构建Docker镜像

```bash
cd BacterialAntigenFinder-2/docker

# 构建完整镜像（含4个conda环境，首次约15-30分钟）
docker compose build
```

> 镜像基于 `continuumio/miniconda3:latest`，内置4个独立conda环境，PyTorch 2.5.1+cu121（兼容CUDA 12.1）。

#### Step 5: 运行示例

**方式A：Docker Compose（推荐）**

```bash
cd docker

# 运行铜绿假单胞菌示例（使用 example_data 中的数据）
docker compose run --rm antigen-finder-test
```

**方式B：docker run 手动运行**

```bash
# 创建输出目录
mkdir -p results/pa

# GPU运行（需要 nvidia-container-toolkit）
docker run --rm --gpus all \
    -v $(pwd)/example_data:/app/data:ro \
    -v $(pwd)/results/pa:/app/results \
    -v $(pwd)/../models:/app/models:ro \
    bacterial-antigen-finder:latest \
    --fasta /app/data/pseudomonas_aeruginosa_antigens.fasta \
    --pdb_dir /app/data/structures/ \
    --organism_type gram- \
    --species Pseudomonas_aeruginosa \
    --models bepipred,discotope,graphbepi,epigraph \
    --output_dir /app/results/ \
    --output_format both
```

> **重要**：BepiPred 模型目录需要**可写**（不加 `:ro`），因为需要缓存 ESM-2 编码。如需写入，将 models 挂载改为 `-v $(pwd)/../models:/app/models`（去掉 `:ro`）。

**仅序列模式（无PDB文件，仅运行BepiPred）**

```bash
docker run --rm --gpus all \
    -v $(pwd)/example_data:/app/data:ro \
    -v $(pwd)/results/pa:/app/results \
    -v $(pwd)/../models:/app/models \
    bacterial-antigen-finder:latest \
    --fasta /app/data/pseudomonas_aeruginosa_antigens.fasta \
    --pdb_dir /app/data/structures/ \
    --organism_type gram- \
    --species Pseudomonas_aeruginosa \
    --models bepipred \
    --output_dir /app/results/ \
    --output_format both
```

#### Step 6: 查看结果

```bash
# 查看输出文件
ls results/pa/
# antigen_candidates.csv  epitope_candidates.csv  analysis_report.html  ...

# 查看候选抗原列表
head -6 results/pa/antigen_candidates.csv

# 查看表位候选列表
head -6 results/pa/epitope_candidates.csv

# 在浏览器中打开HTML报告
# 将 results/pa/analysis_report.html 下载到本地后用浏览器打开
```

#### CPU模式（无GPU）

若无GPU，去掉 `--gpus all`，并添加 `--cpu_only` 参数：

```bash
docker run --rm \
    -v $(pwd)/example_data:/app/data:ro \
    -v $(pwd)/results/pa:/app/results \
    -v $(pwd)/../models:/app/models \
    bacterial-antigen-finder:latest \
    --fasta /app/data/pseudomonas_aeruginosa_antigens.fasta \
    --pdb_dir /app/data/structures/ \
    --species Pseudomonas_aeruginosa \
    --models bepipred \
    --cpu_only \
    --output_dir /app/results/
```

> **注意**：CPU模式下 GraphBepi/EpiGraph 可能因缺少CUDA而无法运行，建议仅使用 `--models bepipred`。

#### GPU挂载方式（NVML故障场景）

如果 `--gpus all` 因 GPU 硬件故障报 NVML 初始化失败，可改用 `--device` 方式绕过：

```bash
# 查询宿主机 NVIDIA 驱动库版本
ls /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.*
ls /lib/x86_64-linux-gnu/libcuda.so.*

# 使用 --device 方式（将 <VERSION> 替换为实际版本号）
docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    --device=/dev/nvidiactl:/dev/nvidiactl \
    --device=/dev/nvidia-uvm:/dev/nvidia-uvm \
    --device=/dev/nvidia-uvm-tools:/dev/nvidia-uvm-tools \
    -v /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.<VERSION>:/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1:ro \
    -v /lib/x86_64-linux-gnu/libcuda.so.<VERSION>:/lib/x86_64-linux-gnu/libcuda.so.1:ro \
    -v $(pwd)/example_data:/app/data:ro \
    -v $(pwd)/results/pa:/app/results \
    -v $(pwd)/../models:/app/models \
    bacterial-antigen-finder:latest \
    --fasta /app/data/pseudomonas_aeruginosa_antigens.fasta \
    --pdb_dir /app/data/structures/ \
    --species Pseudomonas_aeruginosa \
    --output_dir /app/results/
```

### 依赖环境

- Python 3.9+
- PyYAML, pandas, numpy, biopython, plotly
- 各预测模型独立Conda环境（见 `envs/` 目录，Docker镜像已内置）

## 使用方式

### 基本用法

```bash
# 针对铜绿假单胞菌运行
python main.py \
    --fasta example_data/pseudomonas_aeruginosa_antigens.fasta \
    --pdb_dir structures/ \
    --species Pseudomonas_aeruginosa \
    --output_dir results/pa/

# 针对肺炎克雷伯菌运行
python main.py \
    --fasta example_data/klebsiella_pneumoniae_antigens.fasta \
    --pdb_dir structures/ \
    --species Klebsiella_pneumoniae \
    --output_dir results/kp/
```

### 完整参数

```bash
python main.py \
    --fasta input.fasta \
    --pdb_dir structures/ \
    --organism_type gram- \
    --species Pseudomonas_aeruginosa \
    --models bepipred,discotope,graphbepi,epigraph \
    --consensus_threshold 0.5 \
    --min_votes 2 \
    --top_candidates 50 \
    --output_dir results/ \
    --output_format both \
    --workers 4
```

### 参数说明

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--fasta, -f` | str | 必需 | FASTA格式的蛋白质序列文件 |
| `--pdb_dir, -p` | str | 必需 | PDB结构文件目录 |
| `--organism_type, -t` | str | gram- | 生物类型 (gram+, gram-, virus) |
| `--species, -s` | str | None | 目标细菌种类 (如 Pseudomonas_aeruginosa) |
| `--models, -m` | str | 全部 | 启用的模型，逗号分隔 |
| `--consensus_threshold` | float | 0.5 | 共识评分阈值（归一化后分数） |
| `--min_votes` | int | 2 | 最小投票数 |
| `--consensus_method` | str | weighted_avg | 共识方法 (weighted_avg/majority_vote/max_score/mean) |
| `--score_normalize` | str | rank | 跨模型分数归一化 (rank/none) |
| `--top_candidates` | int | 50 | 输出前N个候选 |
| `--output_dir, -o` | str | ./results | 输出目录 |
| `--output_format` | str | both | 输出格式 (csv, html, both) |
| `--cpu_only` | flag | False | 仅使用CPU |
| `--workers` | int | 4 | 并行工作线程数 |
| `--config, -c` | str | None | 自定义配置文件路径 |

### Docker运行

Docker部署和运行的完整指南请参考上方 [方式二：Docker安装](#方式二docker安装推荐用于linux服务器) 章节，包含镜像构建、模型准备、GPU/CPU运行、示例复现等完整步骤。

### 端到端测试

```bash
# 本地Conda环境测试
bash scripts/run_test.sh --species all

# Docker环境测试
bash scripts/run_test_docker.sh --species all

# 输出验证
python scripts/validate_output.py results/pa/
```

## 输入格式

### FASTA序列文件

```fasta
>OprF_PSEAE Pseudomonas aeruginosa Outer membrane protein F
MKLKNTLGVVIGSLVAASAMNAFAQGQNSVEIEAFGKRYFTDSVRNMKNADLYGG...
>OmpK36_KLEPN Klebsiella pneumoniae Porin OmpK36
MKVKVLSLLVPALLVAGAANAAEIYNKDGNKLDLYGKIDGLHYFSDDKSVDGDQT...
```

### PDB结构文件

- 支持标准PDB格式
- 文件名应与FASTA中的序列ID对应（如 `OprF_PSEAE.pdb`）
- 支持实验解析结构和AlphaFold预测结构

## 输出格式

### 抗原候选清单 (antigen_candidates.csv)

蛋白质级别候选，每行一个蛋白质：

| 列名 | 说明 |
|------|------|
| rank | 排名 |
| protein_id | 蛋白质ID |
| epitope_sequence | 最佳表位序列 |
| avg_consensus_score | 平均共识分数 |
| antigenicity_score | VaxiJen抗原性分数 |
| protegenicity_score | 保护性分数（VaxiJen抗原性×100，0-100） |
| composite_score | AI综合评分 |
| recommendation | 推荐等级 (HIGH/MEDIUM/LOW) |

### 表位候选清单 (epitope_candidates.csv)

表位级别候选，每行一个独立表位：

| 列名 | 说明 |
|------|------|
| rank | 表位排名 |
| protein_id | 所属蛋白质 |
| epitope_start / epitope_end | 表位位置 |
| epitope_sequence | 表位氨基酸序列 |
| epitope_length | 表位长度 |
| composite_score | AI综合评分 |
| avg_consensus_score | 共识分数 |
| antigenicity_score | 抗原性分数 |
| protegenicity_score | 保护性分数（VaxiJen抗原性×100） |
| num_predictors_agree | 一致预测器数量 |
| recommendation | 推荐等级 |

### HTML分析报告

- 交互式可视化图表（推荐等级分布、共识分数分布）
- 蛋白质级别 & 表位级别候选表格
- 各预测器统计信息

### 输出目录结构

```
results/
├── antigen_candidates.csv       # 蛋白质级别候选清单
├── epitope_candidates.csv       # 表位级别候选清单
├── analysis_report.html         # HTML分析报告
├── run_metadata.json            # 运行元数据
└── detailed_results/
    ├── consensus_scores.csv     # 残基级共识评分
    └── immunogenicity_scores.csv # 免疫原性评分
```

## 示例结果

我们提供了三种耐药细菌在 **4 模型集成**（BepiPred-3.0 + DiscoTope-3.0 + GraphBepi + EpiGraph）下的完整运行示例，结果位于 `example_output/` 目录：

- `example_output/ecoli/`
- `example_output/klebsiella_pneumoniae/`
- `example_output/pseudomonas_aeruginosa/`

> 默认启用 `--score_normalize rank`，对不同模型的原始分数进行百分位排名归一化后再计算共识分数。启用全部4个模型的指南请参见 `docker/` 目录下的运行示例。

### 运行结果汇总

| 物种 | 中文名 | 输入序列数 | 分析蛋白数 | 共识表位残基 | 抗原候选数 | 表位候选数 | 运行耗时 |
|------|--------|------------|------------|--------------|------------|------------|----------|
| *Escherichia coli* | 大肠杆菌 | 5 | 25 | 1130 | 18 | 50 | ~1432s |
| *Klebsiella pneumoniae* | 肺炎克雷伯菌 | 10 | 25 | 1209 | 19 | 50 | ~1455s |
| *Pseudomonas aeruginosa* | 铜绿假单胞菌 | 10 | 25 | 1159 | 19 | 50 | ~1474s |

> 注：由于 `example_data/structures/` 目录下包含全部三种菌的示例结构，当 `--pdb_dir` 指向该根目录时，Pipeline 会对所有可匹配的结构进行预测，因此单个物种的运行结果中可能包含其它物种的蛋白质。

### Top 5 候选抗原表位（大肠杆菌示例）

| 排名 | 蛋白质ID | 表位位置 | 表位序列 | 长度 | 平均共识评分 | 综合评分 |
|------|----------|----------|----------|------|--------------|----------|
| 1 | OmpA_ECOLI | 312-322 | ASPASNVALLYSGLNARGALAALALEUILEASP | 11 | 0.9331 | 0.7001 |
| 2 | OmpA_ECOLI | 130-136 | ASNVALTYRGLYLYSASNHIS | 7 | 0.9489 | 0.6956 |
| 3 | OmpC_ECOLI | 307-313 | GLYARGGLYTYRASPASPGLU | 7 | 0.8279 | 0.6637 |
| 4 | FimH_ECOLI | 68-74 | ASPTYRPROGLUTHRILETHR | 7 | 0.9016 | 0.6576 |
| 5 | OmpC_ECOLI | 228-232 | GLYARGGLYTYRASP | 5 | 0.8044 | 0.6540 |

> 完整结果（含候选抗原清单、候选表位清单、运行元数据）请查看 `example_output/` 目录。

## 项目结构

```
BacterialAntigenFinder/
├── main.py                              # 主入口CLI
├── config/
│   └── default_config.yaml              # 默认配置（含预定义耐药细菌）
├── src/
│   ├── preprocessor/                    # 数据预处理
│   │   ├── fasta_parser.py              # FASTA解析 & 序列验证
│   │   ├── pdb_validator.py             # PDB结构验证
│   │   └── data_mapper.py               # 序列-结构映射
│   ├── predictors/                      # 预测器包装器
│   │   ├── base_predictor.py            # 预测器基类
│   │   ├── bepipred_wrapper.py          # BepiPred-3.0
│   │   ├── discotope_wrapper.py         # DiscoTope-3.0
│   │   ├── graphbepi_wrapper.py         # GraphBepi
│   │   └── epigraph_wrapper.py          # EpiGraph
│   ├── aggregator/                      # 结果聚合
│   │   ├── consensus_scorer.py          # 共识评分引擎
│   │   ├── immunogenicity.py            # 免疫原性评估
│   │   ├── vaxijen_calculator.py        # VaxiJen 2.0抗原性计算
│   │   └── candidate_ranker.py          # 候选排序（蛋白质级+表位级）
│   └── reporter/                        # 报告生成
│       ├── csv_exporter.py              # CSV导出
│       └── html_report.py               # HTML报告
├── example_data/                        # 示例数据
│   ├── pseudomonas_aeruginosa_antigens.fasta   # 铜绿假单胞菌10个抗原
│   ├── klebsiella_pneumoniae_antigens.fasta    # 肺炎克雷伯菌10个抗原
│   └── sample_antigens.fasta                   # 大肠杆菌示例
├── example_output/                      # 4模型集成运行示例结果
│   ├── ecoli/                           # 大肠杆菌
│   ├── klebsiella_pneumoniae/           # 肺炎克雷伯菌
│   └── pseudomonas_aeruginosa/          # 铜绿假单胞菌
├── examples/                            # 完整运行示例
│   └── docker_gpu_example/              # Docker GPU 多物种示例
│       ├── README.md                    # 示例说明文档
│       ├── ecoli/                       # 大肠杆菌运行结果
│       ├── pseudomonas_aeruginosa/      # 铜绿假单胞菌运行结果
│       └── klebsiella_pneumoniae/       # 肺炎克雷伯菌运行结果
├── docker/                              # Docker部署
│   ├── Dockerfile                       # 完整构建
│   ├── Dockerfile.light                 # 轻量构建
│   ├── docker-compose.yaml              # 编排配置
│   ├── entrypoint.sh                    # 容器入口
│   └── build.sh                         # 构建脚本
├── scripts/                             # 运行脚本
│   ├── setup_envs.sh                    # Conda环境安装
│   ├── run_test.sh                      # 端到端测试
│   ├── run_test_docker.sh              # Docker测试
│   └── validate_output.py              # 输出验证
├── tests/                               # 测试用例
│   ├── test_pipeline.py                 # Pipeline单元测试
│   └── test_dr_resistant_bacteria.py    # 耐药细菌专项测试
├── envs/                                # Conda环境配置
└── docs/                                # 文档
```

## 故障排查

### 1. Docker GPU 相关问题

#### 问题：`docker run --gpus all` 报 NVML 初始化失败

**原因**：服务器 GPU 硬件故障（如某块 GPU 损坏）导致 `nvidia-container-cli` 无法初始化 NVML。

**解决方案**：改用 `--device` 方式直接挂载 GPU 设备文件，绕过 NVML：

```bash
docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    --device=/dev/nvidiactl:/dev/nvidiactl \
    --device=/dev/nvidia-uvm:/dev/nvidia-uvm \
    --device=/dev/nvidia-uvm-tools:/dev/nvidia-uvm-tools \
    -v /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.<VERSION>:/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1:ro \
    -v /lib/x86_64-linux-gnu/libcuda.so.<VERSION>:/lib/x86_64-linux-gnu/libcuda.so.1:ro \
    ... \
    bacterial-antigen-finder:latest
```

#### 问题：PyTorch 报 CUDA 版本不匹配

**原因**：pip 默认安装的 PyTorch 可能使用更高版本的 CUDA（如 12.8），与宿主机驱动不兼容。

**解决方案**：构建镜像时指定 CUDA 12.1 版本的 PyTorch（Dockerfile 已配置）：

```dockerfile
RUN /opt/conda/envs/bepipred_env/bin/pip install --no-cache-dir \
    torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121
```

### 2. BepiPred 相关问题

#### 问题：BepiPred 运行报 `File cannot be opened` 错误

**原因**：BepiPred 需要将 ESM-2 编码缓存到模型目录的 `esm_encodings/` 子目录，若模型目录以只读方式挂载（`:ro`）则无法写入。

**解决方案**：挂载模型目录时**不要**加 `:ro`：

```bash
# 错误 ❌
-v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main:ro

# 正确 ✅
-v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main
```

#### 问题：BepiPred 输出 CSV 列名无法识别

**原因**：BepiPred-3.0 的 `raw_output.csv` 列名为 `Accession, Residue, BepiPred-3.0 score, BepiPred-3.0 linear epitope score`，与早期版本的解析器不兼容。

**解决方案**：已修复，解析器现在自动检测列名。若仍遇到问题，请检查 `src/predictors/bepipred_wrapper.py` 的 `_parse_output` 方法。

### 3. 共识评分相关问题

#### 问题：单模型运行时无共识表位输出

**原因**：默认 `min_votes=2`，单模型运行时 `vote_count` 最大为1，无法满足投票阈值。

**解决方案**：已实现自适应调整机制，`min_votes` 会自动调整为 `min(configured_min_votes, actual_num_predictors)`。若仍遇到问题，可手动指定 `--min_votes 1`。

### 4. 其他问题

#### 问题：`KeyError: 'protein_id'` 当共识评分为空

**解决方案**：已修复，`CandidateRanker` 现在会检查空 DataFrame 并跳过排序。

#### 问题：`_normalize_id` 将 `protein_A.pdb` 识别为 `protein_a` 而非 `protein`

**解决方案**：已修复，现在先移除 `.pdb` 后缀再检查单字母链标识符。

## 引用

如果您使用了本工具，请引用相关模型的论文：

- BepiPred-3.0: Clifford et al., *Protein Science* (2022)
- DiscoTope-3.0: Høie et al., *Frontiers in Immunology* (2024)
- GraphBepi: Zeng et al., *bioRxiv* (2022)
- EpiGraph: Choi & Kim (2023)
- VaxiJen: Doytchinova & Flower, *BMC Bioinformatics* (2007)

## 许可证

本项目采用 MIT 许可证。各子模型请遵循其原始许可证要求。

## 联系方式

如有问题或建议，请提交 [Issue](https://github.com/RuijinHospitalVNAR/AntigenFinder/issues)。
