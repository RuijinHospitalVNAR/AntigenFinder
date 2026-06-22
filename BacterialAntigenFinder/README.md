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

1. **多模型集成预测** — 整合5种B细胞表位预测模型
2. **共识评分** — 加权平均 + 多数投票，降低单一模型偏差
3. **免疫原性评估** — Vaxign-ML保护性抗原预测 + VaxiJen 2.0抗原性计算
4. **双粒度排序** — 蛋白质级别 & 表位级别候选排序
5. **AI综合评分** — 多维度加权评分，输出推荐等级（HIGH/MEDIUM/LOW）
6. **Docker部署** — 完整容器化方案，支持第三方复现

## 实现方式

### 集成模型

| 模型 | 类型 | 输入 | 核心算法 | 权重 |
|------|------|------|----------|------|
| BepiPred-3.0 | 序列基 | FASTA | ESM-2蛋白质语言模型 | 0.25 |
| DiscoTope-3.0 | 结构基 | PDB | ESM-IF1 + 结构特征 | 0.25 |
| GraphBepi | 混合 | PDB | 图神经网络(GNN) | 0.20 |
| EpiGraph | 结构基 | PDB | 图注意力网络(GAT) | 0.20 |
| Vaxign-ML | 序列基 | FASTA | 随机森林 + PSORTb | 0.10 |

### 评分算法

**共识评分**：对每个残基，按模型权重计算加权平均分，结合多数投票（≥2个模型同意）确定表位。

**免疫原性评估**：
- 保护性抗原分数（Vaxign-ML，0-100）
- 抗原性分数（VaxiJen 2.0，基于氨基酸z-scale自相关系数）
- 亚细胞定位权重（外膜/分泌蛋白优先）
- 跨膜区域降权

**综合评分公式**：

蛋白质级别：`0.35×共识分 + 0.25×保护性分 + 0.20×残基数量 + 0.20×区域数量`

表位级别：`0.30×共识分 + 0.25×保护性分 + 0.20×抗原性分 + 0.15×表位长度 + 0.10×模型一致性`

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
│  ├── EpiGraph       ──→  GAT表位预测                            │
│  └── Vaxign-ML      ──→  保护性抗原预测                         │
│       │                                                          │
│       ▼                                                          │
│  Step 3: 共识评分                                                │
│  ├── 加权平均 (权重: 0.25/0.25/0.20/0.20/0.10)                 │
│  ├── 多数投票 (≥2 models agree)                                 │
│  └── 共识表位识别                                                │
│       │                                                          │
│       ▼                                                          │
│  Step 4: 免疫原性评估                                            │
│  ├── Vaxign-ML 保护性抗原评分                                    │
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
cd AntigenFinder/BacterialAntigenFinder

# 2. 运行环境安装脚本
chmod +x scripts/setup_envs.sh
./scripts/setup_envs.sh

# 3. 激活主控环境
conda activate master_env
```

### 方式二：Docker安装（推荐用于GPU服务器）

#### 镜像构建

镜像基于 `continuumio/miniconda3:latest`，内置5个独立 conda 环境，PyTorch 2.5.1+cu121（兼容 CUDA 12.1）。

```bash
# 完整构建（含5个conda环境，约30GB）
cd BacterialAntigenFinder
docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .

# 轻量构建（不含模型权重，需挂载模型目录）
docker build -t bacterial-antigen-finder:light -f docker/Dockerfile.light .
```

#### GPU 挂载方式

本镜像支持两种 GPU 挂载方式：

**方式A：`--device` 直接挂载（推荐，兼容NVML故障场景）**

```bash
docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    --device=/dev/nvidiactl:/dev/nvidiactl \
    --device=/dev/nvidia-uvm:/dev/nvidia-uvm \
    --device=/dev/nvidia-uvm-tools:/dev/nvidia-uvm-tools \
    -v /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.570.86.10:/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1:ro \
    -v /lib/x86_64-linux-gnu/libcuda.so.570.86.10:/lib/x86_64-linux-gnu/libcuda.so.1:ro \
    ... \
    bacterial-antigen-finder:latest
```

**方式B：`--gpus all`（需要 nvidia-container-toolkit 正常工作）**

```bash
docker run --rm --gpus all ... bacterial-antigen-finder:latest
```

> **提示**：若服务器存在 GPU 硬件故障导致 NVML 初始化失败，请使用方式A。`.so` 文件版本号需根据宿主机实际 NVIDIA 驱动版本调整。

#### Docker Compose 编排

```bash
# 使用 GPU 0 运行
GPU_ID=0 docker compose up antigen-finder

# 运行测试服务
docker compose --profile test up antigen-finder-test
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
    --models bepipred,discotope,graphbepi,epigraph,vaxignml \
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
| `--consensus_threshold` | float | 0.5 | 共识评分阈值 |
| `--min_votes` | int | 2 | 最小投票数 |
| `--top_candidates` | int | 50 | 输出前N个候选 |
| `--output_dir, -o` | str | ./results | 输出目录 |
| `--output_format` | str | both | 输出格式 (csv, html, both) |
| `--cpu_only` | flag | False | 仅使用CPU |
| `--workers` | int | 4 | 并行工作线程数 |
| `--config, -c` | str | None | 自定义配置文件路径 |

### Docker运行（GPU加速）

以下示例展示如何使用 Docker GPU 运行铜绿假单胞菌抗原筛选。完整的多物种示例请参考 [examples/docker_gpu_example/](examples/docker_gpu_example/)。

```bash
# 创建输出目录
mkdir -p results/pa

# GPU Docker 运行（--device方式，兼容NVML故障）
docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    --device=/dev/nvidiactl:/dev/nvidiactl \
    --device=/dev/nvidia-uvm:/dev/nvidia-uvm \
    --device=/dev/nvidia-uvm-tools:/dev/nvidia-uvm-tools \
    -v /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.570.86.10:/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1:ro \
    -v /lib/x86_64-linux-gnu/libcuda.so.570.86.10:/lib/x86_64-linux-gnu/libcuda.so.1:ro \
    -v $(pwd)/example_data:/app/data:ro \
    -v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main \
    -v $(pwd)/results/pa:/app/results \
    bacterial-antigen-finder:latest \
    --fasta /app/data/pseudomonas_aeruginosa_antigens.fasta \
    --pdb_dir /app/data/ \
    --organism_type gram- \
    --species Pseudomonas_aeruginosa \
    --models bepipred \
    --output_dir /app/results/ \
    --output_format both \
    --continue_on_error
```

> **重要**：BepiPred 模型目录必须以**读写**方式挂载（不加 `:ro`），因为 BepiPred 需要缓存 ESM-2 编码到 `esm_encodings/` 子目录。

> **提示**：`.so` 文件版本号（如 `570.86.10`）需根据宿主机实际 NVIDIA 驱动版本调整，可通过 `ls /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.*` 查询。

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
| protegenicity_score | 保护性抗原分数 |
| antigenicity_score | VaxiJen抗原性分数 |
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
| protegenicity_score | 保护性分数 |
| antigenicity_score | 抗原性分数 |
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

我们提供了三种耐药细菌的 Docker GPU 完整运行示例，包含候选抗原列表及 AI 综合评分。详见 [examples/docker_gpu_example/](examples/docker_gpu_example/)。

### 运行结果汇总

| 物种 | 中文名 | 输入序列数 | 抗原候选数 | 表位候选数 | 运行耗时 |
|------|--------|------------|------------|------------|----------|
| *Escherichia coli* | 大肠杆菌 | 5 | 5 | 30 | 74.87s |
| *Pseudomonas aeruginosa* | 铜绿假单胞菌 | 10 | 7 | 47 | 90.70s |
| *Klebsiella pneumoniae* | 肺炎克雷伯菌 | 10 | 8 | 47 | 89.34s |

### Top 5 候选抗原表位（铜绿假单胞菌示例）

| 排名 | 蛋白质ID | 表位位置 | 表位序列 | 长度 | 共识评分 | 综合评分 |
|------|----------|----------|----------|------|----------|----------|
| 1 | OprD_PSEAE | 368-389 | DGTKMSDNNVGYKNYGYGEDGK | 22 | 0.5645 | 0.3263 |
| 2 | PcrV_PSEAE | 214-236 | SPKQSGELKGLSDEYPFEKDNNP | 23 | 0.5553 | 0.3254 |
| 3 | PcrV_PSEAE | 165-188 | DAGGIDLVDPTLYGYAVGDPRWKD | 24 | 0.5452 | 0.3242 |
| 4 | OprF_PSEAE | 120-134 | NITNINSDSQGRQQM | 15 | 0.5754 | 0.3137 |
| 5 | OprF_PSEAE | 188-200 | KAAPAPEPVADVC | 13 | 0.5454 | 0.2989 |

### Top 5 候选抗原表位（肺炎克雷伯菌示例）

| 排名 | 蛋白质ID | 表位位置 | 表位序列 | 长度 | 共识评分 | 综合评分 |
|------|----------|----------|----------|------|----------|----------|
| 1 | KpnO_KLEPN | 179-201 | QNAQDINVGTNNRSSDSDVRFDN | 23 | 0.5669 | 0.3289 |
| 2 | OmpK36_KLEPN | 177-197 | NSVSGEGTSPTNNGRGALKQN | 21 | 0.5539 | 0.3212 |
| 3 | OmpA_KLEPN | 31-48 | GFYGNGFQNNNGPTRNDQ | 18 | 0.5587 | 0.3162 |
| 4 | OmpA_KLEPN | 192-210 | EDAAPVVAPAPAPAPEVAT | 19 | 0.5493 | 0.3156 |
| 5 | OmpA_KLEPN | 306-320 | GNTCDNVKARAALID | 15 | 0.5743 | 0.3134 |

> 完整结果（含 HTML 报告、残基级评分）请查看 [examples/docker_gpu_example/](examples/docker_gpu_example/) 目录。

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
│   │   ├── epigraph_wrapper.py          # EpiGraph
│   │   └── vaxignml_wrapper.py          # Vaxign-ML
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
- Vaxign-ML: Ong et al. (2020)
- VaxiJen: Doytchinova & Flower, *BMC Bioinformatics* (2007)

## 许可证

本项目采用 MIT 许可证。各子模型请遵循其原始许可证要求。

## 联系方式

如有问题或建议，请提交 [Issue](https://github.com/RuijinHospitalVNAR/AntigenFinder/issues)。
