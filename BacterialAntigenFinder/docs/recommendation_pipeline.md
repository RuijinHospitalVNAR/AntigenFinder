# 综合推荐流程详解

本文档详细描述BacterialAntigenFinder如何将Vaxign-ML预测结果单独分析，并与共识评分结合，最终输出推荐序列的完整流程。

---

## 整体流程概览

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           综合推荐流程                                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
        ▼                           ▼                           ▼
┌───────────────────┐    ┌───────────────────┐    ┌───────────────────┐
│   表位预测模型     │    │     Vaxign-ML     │    │   蛋白质序列      │
│  (BepiPred等4个)  │    │   保护性抗原预测   │    │   (输入FASTA)    │
└─────────┬─────────┘    └─────────┬─────────┘    └─────────┬─────────┘
          │                        │                        │
          ▼                        ▼                        ▼
┌───────────────────┐    ┌───────────────────┐    ┌───────────────────┐
│  STEP 1: 共识评分  │    │  STEP 2: 免疫原性  │    │   抗原性计算      │
│  (残基级别)        │    │  单独分析          │    │   (序列组成)      │
│                   │    │  (蛋白质级别)      │    │                   │
│  • 加权平均        │    │  • 保护性分数      │    │                   │
│  • 多数投票        │    │  • 亚细胞定位      │    │                   │
│  • 共识表位判定    │    │  • 跨膜区域        │    │                   │
└─────────┬─────────┘    └─────────┬─────────┘    └─────────┬─────────┘
          │                        │                        │
          └────────────────────────┼────────────────────────┘
                                   │
                                   ▼
                    ┌───────────────────────────┐
                    │  STEP 3: 候选排序整合      │
                    │                           │
                    │  • 表位区域提取            │
                    │  • 综合评分计算            │
                    │  • 推荐等级判定            │
                    │  • 最终排序                │
                    └─────────────┬─────────────┘
                                  │
                                  ▼
                    ┌───────────────────────────┐
                    │  STEP 4: 输出推荐序列      │
                    │                           │
                    │  • 候选抗原清单            │
                    │  • 表位序列                │
                    │  • 推荐等级                │
                    └───────────────────────────┘
```

---

## STEP 1: 共识评分 (残基级别)

### 1.1 输入数据

来自4个表位预测模型的残基级预测结果：

| 模型 | 输出 | 权重 |
|------|------|------|
| BepiPred-3.0 | 每个残基的表位概率 (0-1) | 0.25 |
| DiscoTope-3.0 | 每个残基的表位倾向分数 | 0.25 |
| GraphBepi | 每个残基的表位分数 | 0.20 |
| EpiGraph | 每个残基的表位分数 | 0.20 |

### 1.2 共识计算过程

对于每个蛋白质的每个残基：

```python
# 1. 收集各模型预测
residue_predictions = {
    'bepipred': (score=0.82, is_epitope=True),
    'discotope': (score=0.75, is_epitope=True),
    'graphbepi': (score=0.68, is_epitope=True),
    'epigraph': (score=0.45, is_epitope=False)
}

# 2. 计算加权平均
consensus_score = (0.82×0.25 + 0.75×0.25 + 0.68×0.20 + 0.45×0.20) / 
                  (0.25 + 0.25 + 0.20 + 0.20)
               = 0.6856

# 3. 统计投票数
vote_count = 3  # bepipred, discotope, graphbepi投票为表位

# 4. 共识判定
is_consensus_epitope = (consensus_score >= 0.5) AND (vote_count >= 2)
                     = True
```

### 1.3 输出结果

生成 `consensus_scores.csv`：

| protein_id | residue_id | residue_name | consensus_score | vote_count | is_consensus_epitope |
|------------|------------|--------------|-----------------|------------|----------------------|
| OmpA | 45 | K | 0.6856 | 3 | True |
| OmpA | 46 | G | 0.7123 | 4 | True |
| ... | ... | ... | ... | ... | ... |

---

## STEP 2: Vaxign-ML单独分析 (蛋白质级别)

### 2.1 Vaxign-ML预测输出

Vaxign-ML提供蛋白质级别的综合评估：

```
Protein_ID: OmpA_ECOLI
├── Protegenicity Score: 85.2 (百分位)
├── Subcellular Location: OuterMembrane
├── Has Signal Peptide: True
├── Transmembrane Regions: 1
└── Is Protective Antigen: True (score >= 50)
```

### 2.2 免疫原性评估计算

#### 2.2.1 提取Vaxign-ML特征

```python
features = {
    'protegenicity_score': 85.2,
    'subcellular_location': 'OuterMembrane',
    'has_signal_peptide': True,
    'transmembrane_regions': 1
}
```

#### 2.2.2 计算抗原性分数

基于氨基酸组成（类似VaxiJen方法）：

```python
# 氨基酸抗原性系数
antigenicity_coefs = {
    'A': 1.064, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091,
    'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 1.250,
    'M': 0.826, 'N': 0.776, 'P': 1.064, 'Q': 1.015, 'R': 0.873,
    'S': 0.883, 'T': 0.909, 'V': 1.383, 'W': 0.893, 'Y': 1.161
}

# 计算序列的抗原性分数
antigenicity_score = Σ(coef_aa) / sequence_length
# 例: 对于OmpA序列 → antigenicity_score = 1.05
```

#### 2.2.3 计算亚细胞定位权重

表面暴露的蛋白质获得更高权重：

```python
location_weights = {
    'OuterMembrane': 1.5,   # 外膜蛋白 - 最优
    'Fimbrial': 1.5,        # 菌毛蛋白
    'Surface': 1.5,         # 表面蛋白
    'Extracellular': 1.4,   # 胞外蛋白
    'Secreted': 1.4,        # 分泌蛋白
    'Periplasmic': 1.3,     # 周质蛋白
    'CellWall': 1.3,        # 细胞壁蛋白
    'InnerMembrane': 1.0,   # 内膜蛋白
    'Cytoplasmic': 0.8,     # 胞质蛋白 - 降权
    'Unknown': 1.0
}

# OmpA的定位权重
location_weight = 1.5  # OuterMembrane
```

#### 2.2.4 计算跨膜惩罚

过多跨膜区域的蛋白难以表达纯化：

```python
tm_penalty = 0.8 if transmembrane_regions > 2 else 1.0

# OmpA只有1个跨膜区
tm_penalty = 1.0
```

#### 2.2.5 计算综合免疫原性评分

```python
composite_score = protegenicity × location_weight × antigenicity × tm_penalty

# OmpA的计算
composite_score = 85.2 × 1.5 × 1.05 × 1.0 = 134.19
```

### 2.3 推荐等级判定

```python
# 判定逻辑
if is_protective_antigen AND location in PREFERRED_LOCATIONS AND antigenicity >= 0.4:
    recommendation = 'HIGH'
elif is_protective_antigen OR location in PREFERRED_LOCATIONS:
    recommendation = 'MEDIUM'
else:
    recommendation = 'LOW'

# OmpA的判定
# is_protective_antigen = True (85.2 >= 50)
# location = 'OuterMembrane' ∈ PREFERRED_LOCATIONS
# antigenicity = 1.05 >= 0.4
# → recommendation = 'HIGH'
```

### 2.4 输出结果

生成 `immunogenicity_scores.csv`：

| protein_id | protegenicity_score | is_protective | subcellular_location | antigenicity_score | recommendation | composite_score |
|------------|---------------------|---------------|----------------------|--------------------| ---------------|-----------------|
| OmpA | 85.2 | True | OuterMembrane | 1.05 | HIGH | 134.19 |
| FimH | 78.5 | True | Fimbrial | 0.98 | HIGH | 115.40 |
| CysK | 32.1 | False | Cytoplasmic | 0.92 | LOW | 23.67 |

---

## STEP 3: 候选排序整合

### 3.1 提取表位区域

从共识评分结果中提取连续的表位区域：

```python
# 输入: 共识评分DataFrame
# 过滤: is_consensus_epitope == True
# 规则: 连续残基合并为区域，最小长度5

# 示例输出
epitope_regions = [
    (45, 60, avg_score=0.72, sequence="KGFNKTGDTGVLQKA"),
    (120, 135, avg_score=0.68, sequence="NNPVTGQGTANVYVN"),
    (200, 215, avg_score=0.65, sequence="AQGVQLTAKLGYPITD")
]
```

### 3.2 综合评分计算

整合共识评分（表位质量）和免疫原性评估（蛋白质质量）：

```python
# 各分量归一化
score_component = (avg_consensus_score + max_consensus_score) / 2
residue_component = log(1 + total_epitope_residues) / 5
protegen_component = protegenicity_score / 100
region_component = min(num_epitope_regions / 5, 1.0)

# 加权综合
final_composite_score = (
    0.35 × score_component +      # 表位质量权重
    0.25 × protegen_component +   # 保护性分数权重
    0.20 × residue_component +    # 表位残基数量权重
    0.20 × region_component       # 表位区域数量权重
)
```

### 3.3 计算示例

```
蛋白质: OmpA_ECOLI
├── avg_consensus_score: 0.72
├── max_consensus_score: 0.85
├── total_epitope_residues: 42
├── num_epitope_regions: 3
├── protegenicity_score: 85.2
│
├── score_component = (0.72 + 0.85) / 2 = 0.785
├── residue_component = log(1+42) / 5 = 0.748
├── protegen_component = 85.2 / 100 = 0.852
├── region_component = min(3/5, 1.0) = 0.6
│
└── final_score = 0.35×0.785 + 0.25×0.852 + 0.20×0.748 + 0.20×0.6
                = 0.275 + 0.213 + 0.150 + 0.120
                = 0.758
```

### 3.4 最终排序

按 `final_composite_score` 降序排序，生成候选列表。

---

## STEP 4: 输出推荐序列

### 4.1 候选抗原清单 (antigen_candidates.csv)

```csv
rank,protein_id,sequence_length,best_epitope_range,best_epitope_sequence,num_epitope_regions,total_epitope_residues,avg_consensus_score,max_consensus_score,protegenicity_score,subcellular_location,recommendation,composite_score,all_epitope_regions
1,OmpA_ECOLI,346,45-60,KGFNKTGDTGVLQKA,3,42,0.72,0.85,85.2,OuterMembrane,HIGH,0.758,"45-60(0.72); 120-135(0.68); 200-215(0.65)"
2,FimH_ECOLI,300,120-138,NNPVTGQGTANVYVN,2,28,0.68,0.82,78.5,Fimbrial,HIGH,0.712,"120-138(0.68); 245-258(0.65)"
3,LamB_ECOLI,421,85-102,GTSYGYDSQLVQADG,4,55,0.65,0.79,72.3,OuterMembrane,HIGH,0.698,"85-102(0.65); 150-165(0.62); 280-295(0.68); 350-365(0.61)"
```

### 4.2 推荐序列输出格式

```
================================================================================
                    BacterialAntigenFinder - 抗原推荐清单
================================================================================

排名 #1: OmpA_ECOLI
--------------------------------------------------------------------------------
推荐等级: HIGH
综合评分: 0.758
保护性分数: 85.2
亚细胞定位: OuterMembrane
序列长度: 346 aa

表位区域:
  区域1: 残基 45-60 (评分: 0.72)
         序列: KGFNKTGDTGVLQKA
  区域2: 残基 120-135 (评分: 0.68)
         序列: NNPVTGQGTANVYVN
  区域3: 残基 200-215 (评分: 0.65)
         序列: AQGVQLTAKLGYPITD

推荐理由:
  ✓ 保护性抗原 (protegenicity ≥ 50)
  ✓ 优选亚细胞定位 (外膜蛋白)
  ✓ 高抗原性分数
  ✓ 多个高质量表位区域

================================================================================
```

---

## 数据流总结

```
输入序列 (FASTA)
      │
      ├──────────────────────────────────────────────────────┐
      │                                                      │
      ▼                                                      ▼
┌─────────────┐                                    ┌─────────────────┐
│ 4个表位模型  │                                    │    Vaxign-ML    │
│ 残基级预测   │                                    │   蛋白质级预测   │
└──────┬──────┘                                    └────────┬────────┘
       │                                                    │
       ▼                                                    ▼
┌─────────────┐                                    ┌─────────────────┐
│  共识评分    │                                    │   免疫原性评估   │
│ • 加权平均   │                                    │ • 保护性分数     │
│ • 多数投票   │                                    │ • 定位权重      │
│ • 共识表位   │                                    │ • 抗原性分数     │
└──────┬──────┘                                    └────────┬────────┘
       │                                                    │
       │              ┌─────────────────┐                   │
       └──────────────┤   候选排序器    ├───────────────────┘
                      │                 │
                      │ • 表位区域提取  │
                      │ • 综合评分      │
                      │ • 推荐等级      │
                      └────────┬────────┘
                               │
                               ▼
                      ┌─────────────────┐
                      │   推荐序列输出   │
                      │                 │
                      │ • CSV清单       │
                      │ • HTML报告      │
                      │ • 表位序列      │
                      └─────────────────┘
```

---

## 配置参数

在 `config/default_config.yaml` 中可调整：

```yaml
# 共识评分参数
consensus:
  threshold: 0.5          # 共识分数阈值
  min_votes: 2            # 最小投票数
  method: "weighted_avg"  # 评分方法

# 模型权重
models:
  weights:
    bepipred: 0.25
    discotope: 0.25
    graphbepi: 0.20
    epigraph: 0.20
    vaxignml: 0.10       # 用于共识评分

# 免疫原性评估参数
immunogenicity:
  protegenicity_threshold: 50.0   # 保护性抗原阈值
  antigenicity_threshold: 0.4     # 抗原性阈值

# 候选排序参数
ranking:
  min_epitope_length: 5           # 最小表位区域长度
  min_consensus_score: 0.3        # 最小共识分数
  top_candidates: 50              # 输出候选数量

# 综合评分权重
composite_weights:
  score_component: 0.35           # 表位质量
  protegen_component: 0.25        # 保护性分数
  residue_component: 0.20         # 残基数量
  region_component: 0.20          # 区域数量
```

---

## 代码实现位置

| 步骤 | 模块 | 文件 |
|------|------|------|
| 共识评分 | ConsensusScorer | `src/aggregator/consensus_scorer.py` |
| 免疫原性评估 | ImmunogenicityEvaluator | `src/aggregator/immunogenicity.py` |
| 候选排序 | CandidateRanker | `src/aggregator/candidate_ranker.py` |
| 流程整合 | AntigenFinderPipeline | `main.py` |

---

## 关键代码片段

### main.py 中的流程整合

```python
# Step 3: 共识评分
consensus_results = self.consensus_scorer.compute_consensus(predictions)

# Step 4: 免疫原性评估（Vaxign-ML单独分析）
if 'vaxignml' in predictions:
    immunogenicity_scores = self.immunogenicity_evaluator.evaluate(
        predictions['vaxignml'],
        sequences
    )

# Step 5: 候选排序（整合两者）
ranked_candidates = self.candidate_ranker.rank(
    consensus_results,
    immunogenicity_scores,
    sequences
)
```

这确保了Vaxign-ML预测结果经过独立的详细分析后，再与共识评分结合，最终输出高质量的推荐序列。
