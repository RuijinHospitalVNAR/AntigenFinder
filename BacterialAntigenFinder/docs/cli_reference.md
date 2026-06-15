# 命令行参数参考

## 基本用法

```bash
python main.py --fasta INPUT.fasta --pdb_dir STRUCTURES/ --output_dir RESULTS/
```

---

## 完整参数列表

### 输入参数

| 参数 | 短选项 | 类型 | 必需 | 默认值 | 说明 |
|------|--------|------|------|--------|------|
| `--fasta` | `-f` | str | 是 | - | FASTA格式的蛋白质序列文件 |
| `--pdb_dir` | `-p` | str | 是 | - | PDB结构文件目录 |
| `--organism_type` | `-t` | str | 否 | gram- | 生物类型 |

#### organism_type 可选值

| 值 | 说明 |
|----|------|
| `gram-` 或 `g-` | 革兰氏阴性菌 |
| `gram+` 或 `g+` | 革兰氏阳性菌 |
| `virus` 或 `v` | 病毒 |

---

### 模型参数

| 参数 | 短选项 | 类型 | 默认值 | 说明 |
|------|--------|------|--------|------|
| `--models` | `-m` | str | 全部 | 启用的模型，逗号分隔 |
| `--consensus_threshold` | - | float | 0.5 | 共识评分阈值 |
| `--min_votes` | - | int | 2 | 最小投票数 |

#### models 可选值

可以组合使用，用逗号分隔：

```bash
# 使用全部模型
--models bepipred,discotope,graphbepi,epigraph,vaxignml

# 只使用序列基模型
--models bepipred,vaxignml

# 只使用结构基模型
--models discotope,graphbepi,epigraph
```

---

### 输出参数

| 参数 | 短选项 | 类型 | 默认值 | 说明 |
|------|--------|------|--------|------|
| `--output_dir` | `-o` | str | ./results | 输出目录 |
| `--top_candidates` | - | int | 50 | 输出前N个候选 |
| `--output_format` | - | str | both | 输出格式 |

#### output_format 可选值

| 值 | 说明 |
|----|------|
| `csv` | 仅输出CSV文件 |
| `html` | 仅输出HTML报告 |
| `both` | 同时输出CSV和HTML |

---

### 运行参数

| 参数 | 短选项 | 类型 | 默认值 | 说明 |
|------|--------|------|--------|------|
| `--config` | `-c` | str | None | 自定义配置文件路径 |
| `--cpu_only` | - | flag | False | 仅使用CPU |
| `--workers` | - | int | 4 | 并行工作线程数 |
| `--continue_on_error` | - | flag | True | 单个模型失败时继续运行 |
| `--log_level` | - | str | INFO | 日志级别 |
| `--log_file` | - | str | None | 日志文件路径 |

#### log_level 可选值

| 值 | 说明 |
|----|------|
| `DEBUG` | 详细调试信息 |
| `INFO` | 一般信息（推荐） |
| `WARNING` | 警告信息 |
| `ERROR` | 仅错误信息 |

---

## 使用示例

### 示例1: 基本用法

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --output_dir results/
```

### 示例2: 指定生物类型和模型

```bash
python main.py \
    --fasta mrsa_proteins.fasta \
    --pdb_dir mrsa_structures/ \
    --organism_type gram+ \
    --models bepipred,discotope,vaxignml \
    --output_dir mrsa_results/
```

### 示例3: 调整共识参数

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --consensus_threshold 0.6 \
    --min_votes 3 \
    --top_candidates 100 \
    --output_dir results/
```

### 示例4: 仅使用CPU

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --cpu_only \
    --output_dir results/
```

### 示例5: 使用自定义配置

```bash
python main.py \
    --config my_config.yaml \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --output_dir results/
```

### 示例6: 启用详细日志

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --log_level DEBUG \
    --log_file pipeline.log \
    --output_dir results/
```

### 示例7: 快速测试（仅BepiPred）

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --models bepipred \
    --output_format csv \
    --output_dir quick_test/
```

---

## 配置文件格式

除了命令行参数，也可以使用YAML配置文件：

```yaml
# my_config.yaml

# 模型配置
models:
  enabled:
    - bepipred
    - discotope
    - graphbepi
    - epigraph
    - vaxignml
  
  weights:
    bepipred: 0.25
    discotope: 0.25
    graphbepi: 0.20
    epigraph: 0.20
    vaxignml: 0.10
  
  thresholds:
    bepipred: 0.1512
    discotope: 0.90
    graphbepi: 0.1763
    epigraph: 0.1481
    vaxignml: 0.5

# 共识配置
consensus:
  threshold: 0.5
  min_votes: 2
  method: "weighted_avg"

# 输出配置
output:
  top_candidates: 50
  formats:
    - csv
    - html

# 运行配置
runtime:
  use_gpu: true
  max_workers: 4
  timeout: 3600
  continue_on_error: true
```

---

## 帮助信息

查看完整帮助：

```bash
python main.py --help
```

输出：

```
usage: main.py [-h] --fasta FASTA --pdb_dir PDB_DIR 
               [--organism_type {gram+,gram-,g+,g-,virus,v}]
               [--models MODELS] [--consensus_threshold THRESHOLD]
               [--min_votes MIN_VOTES] [--output_dir OUTPUT_DIR]
               [--top_candidates TOP] [--output_format {csv,html,both}]
               [--config CONFIG] [--cpu_only] [--workers WORKERS]
               [--continue_on_error] [--log_level {DEBUG,INFO,WARNING,ERROR}]
               [--log_file LOG_FILE]

BacterialAntigenFinder - 细菌抗原AI智能筛选平台
...
```
