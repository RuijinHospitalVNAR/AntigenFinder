# 快速入门

本指南帮助您在5分钟内完成第一次抗原筛选。

---

## 准备工作

### 1. 激活环境

```bash
conda activate master_env
```

### 2. 准备输入文件

#### FASTA序列文件

创建 `my_antigens.fasta`:

```fasta
>OmpA_ECOLI Outer membrane protein A
MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNG
PTHENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGV
QLTAKLGYPITDDLDVYTRLGGMVWRADTKSNVYGKNHDTGVSPVFAGG

>FimH_ECOLI Type 1 fimbrial adhesin  
MKRVITLFAVLLMGWSVNAWSFACKTANGTAIPIGGGSANVYVNLAPVV
NVGQNLVVDLSTQIFCHNDYPETITDYVTLQRGSAYGGVLSNFSGTVKY
SGSSYPFPTTSETPRVVYNSRTDKPWPVALYLTPVSSAGGVAIKAGSLI
```

#### PDB结构文件目录

```bash
mkdir structures/
# 将对应的PDB文件放入
cp OmpA_ECOLI.pdb structures/
cp FimH_ECOLI.pdb structures/
```

**注意**: PDB文件名应与FASTA中的序列ID匹配。

---

## 运行Pipeline

### 基本命令

```bash
python main.py \
    --fasta my_antigens.fasta \
    --pdb_dir structures/ \
    --organism_type gram- \
    --output_dir my_results/
```

### 预期输出

```
============================================================
BacterialAntigenFinder - 开始抗原筛选流程
版本: 1.0.0
时间: 2024-01-15 10:30:45
============================================================
Step 1: 数据预处理...
  解析到 2 个序列
  验证通过 2 个PDB结构
  成功映射 2 个序列-结构对
Step 2: 运行表位预测模型...
  运行 bepipred...
  bepipred 完成 (15.2s)，300 残基，45 表位
  运行 discotope...
  discotope 完成 (8.5s)，300 残基，52 表位
  ...
Step 3: 计算共识评分...
  共识评分完成，300 个残基
Step 4: 免疫原性评估...
  免疫原性评估完成
Step 5: 候选抗原排序...
  排序完成，Top 2 候选
Step 6: 生成结果报告...
  CSV报告: my_results/antigen_candidates.csv
  HTML报告: my_results/analysis_report.html
============================================================
抗原筛选流程完成!
总耗时: 45.3 秒
结果保存至: my_results/
============================================================
```

---

## 查看结果

### 1. 候选抗原清单

```bash
cat my_results/antigen_candidates.csv
```

或用Excel打开查看。

### 2. 分析报告

用浏览器打开:

```bash
firefox my_results/analysis_report.html
# 或
google-chrome my_results/analysis_report.html
```

### 3. 详细结果

```bash
ls my_results/detailed_results/
# consensus_scores.csv
# immunogenicity_scores.csv
```

---

## 解读结果

### 推荐等级

| 等级 | 含义 | 建议 |
|------|------|------|
| **HIGH** | 高优先级候选 | 优先进行实验验证 |
| **MEDIUM** | 中优先级候选 | 可作为备选 |
| **LOW** | 低优先级候选 | 不推荐 |

### 关键指标

- **Avg_Consensus_Score**: 平均共识分数，≥0.6 为优质
- **Protegenicity_Score**: 保护性分数，≥50 为保护性抗原
- **Subcellular_Location**: 亚细胞定位，优选表面蛋白

---

## 常见场景

### 场景1: 只有序列，没有结构

仅使用序列基模型:

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir empty_dir/ \
    --models bepipred,vaxignml \
    --output_dir results/
```

### 场景2: 快速测试

只用BepiPred快速预览:

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --models bepipred \
    --output_dir quick_test/
```

### 场景3: 更严格的筛选

提高阈值和投票要求:

```bash
python main.py \
    --fasta antigens.fasta \
    --pdb_dir structures/ \
    --consensus_threshold 0.7 \
    --min_votes 3 \
    --output_dir strict_results/
```

### 场景4: 革兰氏阳性菌

```bash
python main.py \
    --fasta gram_positive.fasta \
    --pdb_dir structures/ \
    --organism_type gram+ \
    --output_dir gram_pos_results/
```

---

## 下一步

- [命令行参数详解](cli_reference.md)
- [共识算法原理](consensus_algorithm.md)
- [免疫原性评估](immunogenicity.md)
- [输出格式说明](output_format.md)
