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

---

## Docker 安装（推荐）

此方式已验证在 Linux (Ubuntu 20.04+) 和 Windows (WSL2 + Docker Desktop) 上可成功运行。

### 前置条件

| 需求 | 说明 |
|------|------|
| **操作系统** | Linux（推荐 Ubuntu 20.04+）或 Windows WSL2 + Docker Desktop |
| **Docker** | 20.10+、Docker Compose v2 |
| **GPU（可选）** | NVIDIA GPU + 驱动 + [nvidia-container-toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) |
| **内存** | GPU 模式 ≥ 8GB；CPU 模式 ≥ 16GB（仅 BepiPred + DiscoTope 需 8GB） |
| **磁盘** | ≥ 60GB（镜像约 30GB + 模型文件约 10GB + 示例数据） |

### Step 1: 克隆项目

```bash
git clone https://github.com/RuijinHospitalVNAR/AntigenFinder.git
cd AntigenFinder
```

### Step 2: 准备预测模型

4 个预测模型需要克隆到 `AntigenFinder/` 目录下（与 `BacterialAntigenFinder-2/` 同级）：

```bash
# 在 AntigenFinder/ 目录下
git clone https://github.com/mabu-car/BepiPred-3.0.git BepiPred-3.0-main
git clone https://github.com/mabu-car/DiscoTope-3.0.git DiscoTope-3.0-master
git clone https://github.com/GraphBepi/GraphBepi.git GraphBepi-main
git clone https://github.com/GraphBepi/EpiGraph.git EpiGraph-main
```

> **注意**：BepiPred 模型运行时会在 `esm_encodings/` 目录下缓存 ESM-2 编码结果，因此**不要**对该目录设置只读权限。

### Step 2b: 解压 DiscoTope 模型权重

DiscoTope-3.0 的 XGBoost 模型权重被压缩在 `models.zip` 中，需要手动解压：

```bash
cd DiscoTope-3.0-master
# Linux
unzip models.zip -d models/
# 如果解压后嵌套了 models/models/ 目录，手动移出
if [ -d models/models ]; then
    mv models/models/* models/ && rmdir models/models
fi
cd ..
```

### Step 3: 预下载 ESM 模型文件（重要）

BepiPred、GraphBepi、EpiGraph 都依赖 Meta 的 ESM 蛋白质语言模型权重文件。这些文件从 `dl.fbaipublicfiles.com` 下载，**国内用户可能无法直接访问**，建议提前下载。

```bash
mkdir -p esm_cache && cd esm_cache

# 方式一：直接下载（海外服务器）
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t36_3B_UR50D.pt
wget https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt
wget https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t36_3B_UR50D-contact-regression.pt
# GraphBepi/EpiGraph 需要
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt

# 方式二：HuggingFace 镜像（国内推荐）
wget https://hf-mirror.com/facebook/esm2_t33_650M_UR50D/resolve/main/model.pt \
    -O esm2_t33_650M_UR50D.pt
wget https://hf-mirror.com/facebook/esm2_t36_3B_UR50D/resolve/main/model.pt \
    -O esm2_t36_3B_UR50D.pt
wget https://hf-mirror.com/facebook/esm2_t33_650M_UR50D/resolve/main/contact_regression.pt \
    -O esm2_t33_650M_UR50D-contact-regression.pt
wget https://hf-mirror.com/facebook/esm2_t36_3B_UR50D/resolve/main/contact_regression.pt \
    -O esm2_t36_3B_UR50D-contact-regression.pt
wget https://hf-mirror.com/facebook/esm_if1_gvp4_t16_142M_UR50/resolve/main/model.pt \
    -O esm_if1_gvp4_t16_142M_UR50.pt

cd ..
```

| 文件 | 大小 | 被谁使用 |
|------|------|----------|
| `esm2_t33_650M_UR50D.pt` | ~2.5 GB | BepiPred |
| `esm2_t33_650M_UR50D-contact-regression.pt` | ~4 KB | BepiPred |
| `esm2_t36_3B_UR50D.pt` | ~5.4 GB | GraphBepi |
| `esm2_t36_3B_UR50D-contact-regression.pt` | ~7 KB | GraphBepi |
| `esm_if1_gvp4_t16_142M_UR50.pt` | ~1.6 GB | GraphBepi, EpiGraph |

### Step 4: 目录结构确认

完成上述步骤后，确保如下结构：

```
AntigenFinder/
├── BepiPred-3.0-main/           # BepiPred 模型（需可写）
├── DiscoTope-3.0-master/        # DiscoTope 模型（已解压 models.zip）
│   └── models/                  # 含 100+ 个 XGBoost JSON 文件
├── GraphBepi-main/              # GraphBepi 模型
├── EpiGraph-main/               # EpiGraph 模型
├── esm_cache/                   # ESM 预下载权重（映射到容器内 /root/.cache/torch/hub/checkpoints）
│   ├── esm2_t33_650M_UR50D.pt
│   ├── esm2_t36_3B_UR50D.pt
│   ├── esm_if1_gvp4_t16_142M_UR50.pt
│   └── ...-regression.pt
└── BacterialAntigenFinder-2/    # 项目代码
    ├── docker/
    │   ├── Dockerfile
    │   ├── docker-compose.yaml
    │   └── entrypoint.sh
    ├── example_data/
    ├── envs/
    └── src/
```

### Step 5: 构建 Docker 镜像

```bash
cd BacterialAntigenFinder-2

# Docker Compose 构建（注意：compose 文件在 docker/ 子目录）
docker compose -f docker/docker-compose.yaml build

# 或直接 docker build
docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .
```

> **关于镜像源**：`docker/Dockerfile` 默认配置了清华/交大镜像加速（conda、pip、PyTorch），适合国内网络环境。海外服务器可去掉 `Dockerfile` 中第 27-42 行的镜像源配置以加速构建。
>
> 构建时间：首次约 15-30 分钟（需下载 miniconda3 基础镜像 + 创建 5 个 conda 环境）。后续增量构建利用 Docker 层缓存，仅修改的文件层需要重新构建。

### Step 6: 运行示例

**GPU 模式（推荐，Linux + NVIDIA）：**

```bash
# 大肠杆菌示例数据测试
docker run --rm --gpus 1 \
    -v $(pwd)/example_data:/app/data:ro \
    -v $(pwd)/results/ecoli:/app/results \
    -v $(pwd)/../BepiPred-3.0-main:/app/models/BepiPred-3.0-main \
    -v $(pwd)/../DiscoTope-3.0-master:/app/models/DiscoTope-3.0-master:ro \
    -v $(pwd)/../GraphBepi-main:/app/models/GraphBepi-main:ro \
    -v $(pwd)/../EpiGraph-main:/app/models/EpiGraph-main:ro \
    -v $(pwd)/../esm_cache:/root/.cache/torch/hub/checkpoints:ro \
    bacterial-antigen-finder:latest \
    --fasta /app/data/sample_antigens.fasta \
    --pdb_dir /app/data/structures/ \
    --organism_type gram- \
    --species Escherichia_coli \
    --models bepipred,discotope,graphbepi,epigraph \
    --output_dir /app/results/ \
    --output_format both
```

**CPU 模式（无 GPU，仅 BepiPred + DiscoTope）：**

```bash
docker run --rm \
    -v $(pwd)/example_data:/app/data:ro \
    -v $(pwd)/results/ecoli:/app/results \
    -v $(pwd)/../BepiPred-3.0-main:/app/models/BepiPred-3.0-main \
    -v $(pwd)/../DiscoTope-3.0-master:/app/models/DiscoTope-3.0-master:ro \
    -v $(pwd)/../esm_cache:/root/.cache/torch/hub/checkpoints:ro \
    bacterial-antigen-finder:latest \
    --fasta /app/data/sample_antigens.fasta \
    --pdb_dir /app/data/structures/ \
    --organism_type gram- \
    --species Escherichia_coli \
    --models bepipred,discotope \
    --output_dir /app/results/ \
    --output_format both \
    --cpu_only
```

> **注意**：
> - BepiPred 模型目录**不加 `:ro`**（需写入 ESM 编码缓存）
> - CPU 模式下 GraphBepi/EpiGraph 需要大量内存（3B 参数模型），Docker 内存 ≥ 16GB 时可能勉强运行，但每个结构耗时 10-30 分钟。建议 CPU 模式只使用 `--models bepipred,discotope`
> - 使用 Docker Compose 时，`MODEL_BASE_DIR` 环境变量可指定模型目录的父路径，`ESM_CACHE_DIR` 可指定 ESM 缓存目录

### Step 7: 查看结果

```bash
# 查看输出文件
ls results/ecoli/
# analysis_report.html  antigen_candidates.csv  epitope_candidates.csv  run_metadata.json

# 候选抗原
head -6 results/ecoli/antigen_candidates.csv

# HTML 报告（浏览器打开）
# results/ecoli/analysis_report.html
```

### 故障排查

#### 构建阶段

| 问题 | 原因 | 解决 |
|------|------|------|
| conda `IncompleteRead` / HTTP 错误 | 镜像源不稳定 | Dockerfile 已配清华镜像 + 重试机制。如仍失败，重新构建（利用缓存从断点继续） |
| PyTorch 下载失败 | `download.pytorch.org` 不可达 | Dockerfile 已配交大镜像。海外服务器去掉 Dockerfile 第 68 行的 `-f` 参数 |
| pip 依赖下载超时 | 网络不稳定 | 已配置 `PIP_RETRIES=10`，重试即可 |
| `graphbepi_env` 报 `repo.anaconda.com` 错误 | env yaml 文件中 `defaults` 频道覆盖了清华镜像 | 已从所有 `envs/*.yaml` 移除 `- defaults` 行 |

#### 运行阶段

| 问题 | 原因 | 解决 |
|------|------|------|
| BepiPred 报 `File cannot be opened` | 模型目录以 `:ro` 挂载，无法写入 ESM 编码 | 去掉 `:ro` 标志 |
| DiscoTope 报 `No module named 'discotope3'` | 已修复（`conda run` PYTHONPATH 传递问题） | 确保使用最新代码 |
| DiscoTope 报 `Found 0/100 XGBoost model JSON files` | 未解压 `models.zip` | 执行 Step 2b |
| ESM 下载报 `dl.fbaipublicfiles.com` 连接超时 | 该域名国内不可达 | 预下载模型到 `esm_cache/`（Step 3） |
| GraphBepi/EpiGraph 报 `Killed` | 内存不足 (OOM) | 增加 Docker 内存至 16GB+ 或使用 GPU |

---

## 本地 Conda 环境安装（备选）

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

### 依赖环境

- Python 3.9+
- PyYAML, pandas, numpy, biopython, plotly
- 各预测模型独立Conda环境（见 `envs/` 目录）

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
├── envs/                                # Conda环境配置
└── docs/                                # 文档
```

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
