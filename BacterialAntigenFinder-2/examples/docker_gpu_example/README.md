# Docker GPU 多物种抗原筛选示例

本示例展示如何使用 AntigenFinder Docker 镜像（GPU 加速）对三种耐药细菌进行抗原筛选，获取候选抗原列表及 AI 综合评分。

> **模型说明**：本示例使用 **BepiPred-3.0**（基于 ESM-2 蛋白质语言模型的序列基 B细胞表位预测，权重 0.25）作为最简示例。如需启用全部4个模型（DiscoTope-3.0/GraphBepi/EpiGraph），请提供 PDB 结构文件并参见下方"进阶：启用全部模型"章节。最新4模型集成结果已预置于 `example_output/` 目录。

---

## 示例概述

| 物种 | 中文名 | 耐药类型 | 输入序列数 | 抗原候选数 | 表位候选数 | 共识残基数 | 运行耗时 |
|------|--------|----------|------------|------------|------------|------------|----------|
| *Escherichia coli* | 大肠杆菌 | ESBL | 5 | 5 | 30 | 1446 | 74.87s |
| *Pseudomonas aeruginosa* | 铜绿假单胞菌 | MDR/XDR | 10 | 7 | 47 | 3264 | 90.70s |
| *Klebsiella pneumoniae* | 肺炎克雷伯菌 | ESBL/CRE | 10 | 8 | 47 | 3264 | 89.34s |

**运行环境**：
- Docker 镜像：`bacterial-antigen-finder:latest`（基于 `continuumio/miniconda3:latest`）
- GPU：NVIDIA（通过 `--device` 方式挂载，兼容 NVML 故障场景）
- PyTorch：2.5.1+cu121
- 启用模型：BepiPred-3.0（序列基，ESM-2 蛋白质语言模型）
- 阈值：0.1512（BepiPred 默认），共识阈值 0.5，min_votes 自适应为 1

---

## 环境要求

1. **Docker**：已安装并可正常运行
2. **NVIDIA GPU**：驱动版本 ≥ 535，CUDA 兼容 12.1
3. **NVIDIA 库文件**：宿主机需存在以下文件（版本号根据实际驱动调整）：
   - `/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.570.86.10`
   - `/lib/x86_64-linux-gnu/libcuda.so.570.86.10`
4. **BepiPred 模型**：`BepiPred-3.0-main` 目录（需可读写，用于缓存 ESM encodings）

---

## 完整运行命令

### 1. 大肠杆菌示例

```bash
mkdir -p results/ecoli

docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    --device=/dev/nvidiactl:/dev/nvidiactl \
    --device=/dev/nvidia-uvm:/dev/nvidia-uvm \
    --device=/dev/nvidia-uvm-tools:/dev/nvidia-uvm-tools \
    -v /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.570.86.10:/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1:ro \
    -v /lib/x86_64-linux-gnu/libcuda.so.570.86.10:/lib/x86_64-linux-gnu/libcuda.so.1:ro \
    -v $(pwd)/example_data:/app/data:ro \
    -v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main \
    -v $(pwd)/results/ecoli:/app/results \
    bacterial-antigen-finder:latest \
    --fasta /app/data/sample_antigens.fasta \
    --pdb_dir /app/data/ \
    --organism_type gram- \
    --species Escherichia_coli \
    --models bepipred \
    --output_dir /app/results/ \
    --output_format both \
    --continue_on_error
```

### 2. 铜绿假单胞菌示例

```bash
mkdir -p results/pseudomonas_aeruginosa

docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    --device=/dev/nvidiactl:/dev/nvidiactl \
    --device=/dev/nvidia-uvm:/dev/nvidia-uvm \
    --device=/dev/nvidia-uvm-tools:/dev/nvidia-uvm-tools \
    -v /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.570.86.10:/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1:ro \
    -v /lib/x86_64-linux-gnu/libcuda.so.570.86.10:/lib/x86_64-linux-gnu/libcuda.so.1:ro \
    -v $(pwd)/example_data:/app/data:ro \
    -v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main \
    -v $(pwd)/results/pseudomonas_aeruginosa:/app/results \
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

### 3. 肺炎克雷伯菌示例

```bash
mkdir -p results/klebsiella_pneumoniae

docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    --device=/dev/nvidiactl:/dev/nvidiactl \
    --device=/dev/nvidia-uvm:/dev/nvidia-uvm \
    --device=/dev/nvidia-uvm-tools:/dev/nvidia-uvm-tools \
    -v /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.570.86.10:/usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1:ro \
    -v /lib/x86_64-linux-gnu/libcuda.so.570.86.10:/lib/x86_64-linux-gnu/libcuda.so.1:ro \
    -v $(pwd)/example_data:/app/data:ro \
    -v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main \
    -v $(pwd)/results/klebsiella_pneumoniae:/app/results \
    bacterial-antigen-finder:latest \
    --fasta /app/data/klebsiella_pneumoniae_antigens.fasta \
    --pdb_dir /app/data/ \
    --organism_type gram- \
    --species Klebsiella_pneumoniae \
    --models bepipred \
    --output_dir /app/results/ \
    --output_format both \
    --continue_on_error
```

---

## Top 5 候选抗原表位（按 AI 综合评分排序）

### 大肠杆菌 (*Escherichia coli*)

| 排名 | 蛋白质ID | 表位位置 | 表位序列 | 长度 | 共识评分 | 综合评分 | 推荐等级 |
|------|----------|----------|----------|------|----------|----------|----------|
| 1 | OmpC_ECOLI | 185-223 | RDTARRSNGDGVGGSITYEVTGKNGSQSQPYSGNTNGAI | 39 | 0.5608 | 0.3383 | LOW |
| 2 | OmpC_ECOLI | 84-106 | YNFQGNNSEGADAQTGNKTRLAF | 23 | 0.5963 | 0.3377 | LOW |
| 3 | LamB_ECOLI | 172-200 | SEAGGSSFATYERYLAADGAWRHDAETGT | 29 | 0.5409 | 0.3308 | LOW |
| 4 | OmpC_ECOLI | 251-271 | YNFQGNNSEGADAQTGNKTRL | 21 | 0.5843 | 0.3303 | LOW |
| 5 | OmpA_ECOLI | 194-212 | GEAAPVVAPAPAPAPEVQT | 19 | 0.5518 | 0.3164 | LOW |

### 铜绿假单胞菌 (*Pseudomonas aeruginosa*)

| 排名 | 蛋白质ID | 表位位置 | 表位序列 | 长度 | 共识评分 | 综合评分 | 推荐等级 |
|------|----------|----------|----------|------|----------|----------|----------|
| 1 | OprD_PSEAE | 368-389 | DGTKMSDNNVGYKNYGYGEDGK | 22 | 0.5645 | 0.3263 | LOW |
| 2 | PcrV_PSEAE | 214-236 | SPKQSGELKGLSDEYPFEKDNNP | 23 | 0.5553 | 0.3254 | LOW |
| 3 | PcrV_PSEAE | 165-188 | DAGGIDLVDPTLYGYAVGDPRWKD | 24 | 0.5452 | 0.3242 | LOW |
| 4 | OprF_PSEAE | 120-134 | NITNINSDSQGRQQM | 15 | 0.5754 | 0.3137 | LOW |
| 5 | OprF_PSEAE | 188-200 | KAAPAPEPVADVC | 13 | 0.5454 | 0.2989 | LOW |

### 肺炎克雷伯菌 (*Klebsiella pneumoniae*)

| 排名 | 蛋白质ID | 表位位置 | 表位序列 | 长度 | 共识评分 | 综合评分 | 推荐等级 |
|------|----------|----------|----------|------|----------|----------|----------|
| 1 | KpnO_KLEPN | 179-201 | QNAQDINVGTNNRSSDSDVRFDN | 23 | 0.5669 | 0.3289 | LOW |
| 2 | OmpK36_KLEPN | 177-197 | NSVSGEGTSPTNNGRGALKQN | 21 | 0.5539 | 0.3212 | LOW |
| 3 | OmpA_KLEPN | 31-48 | GFYGNGFQNNNGPTRNDQ | 18 | 0.5587 | 0.3162 | LOW |
| 4 | OmpA_KLEPN | 192-210 | EDAAPVVAPAPAPAPEVAT | 19 | 0.5493 | 0.3156 | LOW |
| 5 | OmpA_KLEPN | 306-320 | GNTCDNVKARAALID | 15 | 0.5743 | 0.3134 | LOW |

---

## 蛋白质级别候选抗原（按综合评分排序）

### 大肠杆菌 — Top 5 抗原候选

| 排名 | 蛋白质ID | 序列长度 | 表位区域数 | 表位残基总数 | 平均共识分 | 最高共识分 | 综合评分 |
|------|----------|----------|------------|--------------|------------|------------|----------|
| 1 | OmpC_ECOLI | 290 | 7 | 169 | 0.5589 | 0.6303 | 0.6135 |
| 2 | OmpA_ECOLI | 346 | 5 | 135 | 0.5460 | 0.6339 | 0.6030 |
| 3 | PhoE_ECOLI | 300 | 7 | 123 | 0.5513 | 0.6066 | 0.5954 |
| 4 | FimH_ECOLI | 300 | 6 | 115 | 0.5360 | 0.6200 | 0.5924 |
| 5 | LamB_ECOLI | 210 | 5 | 92 | 0.5424 | 0.6057 | 0.5822 |

### 铜绿假单胞菌 — Top 5 抗原候选

| 排名 | 蛋白质ID | 序列长度 | 表位区域数 | 表位残基总数 | 平均共识分 | 最高共识分 | 综合评分 |
|------|----------|----------|------------|--------------|------------|------------|----------|
| 1 | ToxA_PSEAE | 638 | 16 | 266 | 0.5265 | 0.6007 | 0.6208 |
| 2 | OprF_PSEAE | 350 | 8 | 130 | 0.5396 | 0.6303 | 0.5997 |
| 3 | OprD_PSEAE | 443 | 8 | 122 | 0.5353 | 0.6060 | 0.5922 |
| 4 | PcrV_PSEAE | 294 | 4 | 129 | 0.5326 | 0.5941 | 0.5519 |
| 5 | LasB_PSEAE | 498 | 4 | 100 | 0.5356 | 0.6286 | 0.5483 |

### 肺炎克雷伯菌 — Top 5 抗原候选

| 排名 | 蛋白质ID | 序列长度 | 表位区域数 | 表位残基总数 | 平均共识分 | 最高共识分 | 综合评分 |
|------|----------|----------|------------|--------------|------------|------------|----------|
| 1 | KpnO_KLEPN | 374 | 8 | 148 | 0.5499 | 0.6203 | 0.6049 |
| 2 | OmpA_KLEPN | 344 | 6 | 140 | 0.5450 | 0.6340 | 0.6043 |
| 3 | OmpK35_KLEPN | 359 | 10 | 142 | 0.5431 | 0.6087 | 0.6001 |
| 4 | OmpK36_KLEPN | 367 | 7 | 136 | 0.5457 | 0.6108 | 0.5992 |
| 5 | FimH_KLEPN | 302 | 7 | 119 | 0.5362 | 0.6218 | 0.5941 |

---

## 输出文件说明

每个物种子目录包含以下文件：

| 文件 | 说明 |
|------|------|
| `antigen_candidates.csv` | 蛋白质级别候选抗原清单（含综合评分、表位区域、推荐等级） |
| `epitope_candidates.csv` | 表位级别候选清单（含表位序列、位置、共识评分、综合评分） |
| `analysis_report.html` | HTML 可视化分析报告（推荐等级分布、共识分数分布、候选表格） |
| `run_metadata.json` | 运行元数据（版本、耗时、输入输出统计） |
| `detailed_results/consensus_scores.csv` | 残基级共识评分（每个残基的各模型评分与投票情况） |

---

## 进阶：启用全部模型

AntigenFinder 设计支持4个预测模型（权重合计 1.0）：

| 模型 | 类型 | 需要PDB | 权重 | 状态 |
|------|------|---------|------|------|
| BepiPred-3.0 | 序列（ESM-2） | 否 | 0.25 | 已启用 |
| DiscoTope-3.0 | 结构（ESM-IF1） | **是** | 0.25 | 需PDB |
| GraphBepi | 图神经网络 | **是** | 0.25 | 需PDB |
| EpiGraph | 图注意力网络 | **是** | 0.25 | 需PDB |

### 启用全部模型的方法

**1. 准备 PDB 结构文件**

将每个蛋白质的 PDB 结构文件放入 `pdb_dir` 目录。可从以下来源获取：
- [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)
- [RCSB PDB](https://www.rcsb.org/)
- 本地 AlphaFold2 / ESMFold 预测

**2. 运行命令**

```bash
docker run --rm \
    --device=/dev/nvidia0:/dev/nvidia0 \
    ... \
    -v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main \
    -v /path/to/DiscoTope-3.0-master:/app/models/DiscoTope-3.0-master \
    -v /path/to/GraphBepi-main:/app/models/GraphBepi-main \
    -v /path/to/EpiGraph-main:/app/models/EpiGraph-main \
    -v $(pwd)/pdb_structures:/app/data:ro \
    -v $(pwd)/results:/app/results \
    bacterial-antigen-finder:latest \
    --fasta /app/data/proteins.fasta \
    --pdb_dir /app/data/ \
    --organism_type gram- \
    --species Pseudomonas_aeruginosa \
    --models bepipred,discotope,graphbepi,epigraph \
    --output_dir /app/results/ \
    --continue_on_error
```



### 多模型效果

启用多模型后，共识评分默认基于百分位排名归一化后的加权平均（BepiPred 0.25 + DiscoTope 0.25 + GraphBepi 0.25 + EpiGraph 0.25），结合多数投票机制，可显著提升表位预测的准确性和可靠性。

---

## 复现步骤

1. **构建 Docker 镜像**（如未构建）：
   ```bash
   cd BacterialAntigenFinder
   docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .
   ```

2. **准备 BepiPred 模型目录**：
   ```bash
   # 确保模型目录可读写（BepiPred 需缓存 ESM encodings）
   chmod -R u+w /path/to/BepiPred-3.0-main
   ```

3. **运行三个物种的 Pipeline**（按上述命令依次执行）

4. **验证输出**：
   ```bash
   # 检查每个物种的候选数量
   for sp in ecoli pseudomonas_aeruginosa klebsiella_pneumoniae; do
       echo "=== $sp ==="
       wc -l results/$sp/antigen_candidates.csv
       wc -l results/$sp/epitope_candidates.csv
   done
   ```

5. **查看 HTML 报告**：
   ```bash
   # 在浏览器中打开
   xdg-open results/pseudomonas_aeruginosa/analysis_report.html
   ```

---

## 注意事项

1. **模型挂载必须可读写**：BepiPred-3.0 需要将 ESM-2 编码缓存到 `esm_encodings/` 子目录，因此挂载时**不要**使用 `:ro`。
2. **GPU 设备挂载**：本示例使用 `--device` 方式而非 `--gpus all`，兼容 NVML 故障的服务器。若 `nvidia-container-toolkit` 正常，也可使用 `--gpus all`。
3. **单模型运行**：本示例仅启用 BepiPred（序列基模型）。`min_votes` 会自适应调整为 1，确保单模型也能产出共识表位。
4. **完整模型运行**：如需启用全部4个模型，需提供对应的 PDB 结构文件，并通过 `--models bepipred,discotope,graphbepi,epigraph` 指定。默认已启用 `--score_normalize rank` 以处理不同模型的分数尺度差异。
5. **`.so` 文件版本**：示例中使用 `570.86.10`，请根据宿主机实际 NVIDIA 驱动版本调整。
