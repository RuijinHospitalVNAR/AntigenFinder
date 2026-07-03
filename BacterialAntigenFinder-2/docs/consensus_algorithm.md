# 共识评分算法详解

## 概述

共识评分模块(`ConsensusScorer`)整合5个预测模型的结果，通过加权平均和多数投票双重机制，计算每个残基的共识表位分数。

## 算法流程

```
┌──────────────────────────────────────────────────────────────┐
│  输入: 5个模型的预测结果                                       │
│  {model_name: PredictionResult}                               │
└──────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│  Step 1: 收集残基预测                                         │
│  按(protein_id, residue_id)分组，收集各模型对该残基的预测      │
│                                                              │
│  {(protein_id, residue_id): {                                │
│      'bepipred': (score=0.8, is_epitope=True),               │
│      'discotope': (score=0.7, is_epitope=True),              │
│      ...                                                     │
│  }}                                                          │
└──────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│  Step 2: 计算加权平均分数                                     │
│                                                              │
│                 Σ(score_i × weight_i)                        │
│  consensus = ─────────────────────────                       │
│                    Σ(weight_i)                               │
│                                                              │
│  默认权重:                                                    │
│    bepipred=0.25, discotope=0.25,                           │
│    graphbepi=0.20, epigraph=0.20, vaxignml=0.10             │
└──────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│  Step 3: 统计投票数                                           │
│                                                              │
│  vote_count = Σ(1 if model predicts epitope else 0)         │
│                                                              │
│  例: 如果bepipred、discotope、graphbepi都预测为表位          │
│      → vote_count = 3                                        │
└──────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│  Step 4: 共识表位判定                                         │
│                                                              │
│  is_consensus_epitope = (consensus_score >= threshold)       │
│                         AND (vote_count >= min_votes)        │
│                                                              │
│  默认: threshold=0.5, min_votes=2                            │
└──────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│  输出: ConsensusResult                                        │
│  protein_id, residue_id, residue_name,                       │
│  consensus_score, vote_count, is_consensus_epitope           │
└──────────────────────────────────────────────────────────────┘
```

## 权重分配原则

### 默认权重

```python
DEFAULT_WEIGHTS = {
    'bepipred': 0.25,    # 序列基，ESM-2语言模型
    'discotope': 0.25,   # 结构基，ESM-IF1逆折叠
    'graphbepi': 0.20,   # GNN，需要结构
    'epigraph': 0.20,    # GAT，需要结构
    'vaxignml': 0.10     # 蛋白质级别，辅助参考
}
```

### 权重设计理由

| 模型 | 权重 | 理由 |
|------|------|------|
| BepiPred-3.0 | 0.25 | 最新的序列基方法，基于ESM-2，在基准测试中表现优异 |
| DiscoTope-3.0 | 0.25 | 最新的结构基方法，利用ESM-IF1逆折叠表示 |
| GraphBepi | 0.20 | 图神经网络方法，捕获空间邻近关系 |
| EpiGraph | 0.20 | 图注意力网络，10折集成提高鲁棒性 |
| Vaxign-ML | 0.10 | 蛋白质级别预测，主要用于免疫原性评估 |

## 计算示例

### 示例数据

假设蛋白质`OmpA`的第45号残基(Ala)被各模型预测如下：

| 模型 | 分数 | 是否表位 | 权重 |
|------|------|----------|------|
| bepipred | 0.82 | True | 0.25 |
| discotope | 0.75 | True | 0.25 |
| graphbepi | 0.68 | True | 0.20 |
| epigraph | 0.45 | False | 0.20 |
| vaxignml | 0.60 | True | 0.10 |

### 计算过程

**1. 加权平均分数**

```
consensus = (0.82×0.25 + 0.75×0.25 + 0.68×0.20 + 0.45×0.20 + 0.60×0.10)
          / (0.25 + 0.25 + 0.20 + 0.20 + 0.10)

         = (0.205 + 0.1875 + 0.136 + 0.09 + 0.06) / 1.0
         = 0.6785
```

**2. 投票统计**

```
vote_count = 4  (bepipred, discotope, graphbepi, vaxignml)
```

**3. 共识判定**

```
consensus_score (0.6785) >= threshold (0.5) → True
vote_count (4) >= min_votes (2) → True

→ is_consensus_epitope = True
```

## 可选评分方法

### 1. 加权平均 (weighted_avg) - 默认

```python
consensus_score = Σ(score_i × weight_i) / Σ(weight_i)
```

### 2. 多数投票 (majority_vote)

```python
consensus_score = vote_count / total_models
```

只看"是否为表位"的二元判断，忽略具体分数。

### 3. 最大分数 (max_score)

```python
consensus_score = max(score_i for all models)
```

取所有模型中的最高分，适合发现潜在表位。

## 表位区域提取

### 连续区域定义

共识表位不是单个残基，而是连续的残基区域。提取规则：

1. 收集所有 `is_consensus_epitope=True` 的残基
2. 按残基ID排序
3. 合并连续残基（残基ID差值=1）
4. 过滤长度 < `min_length`（默认5）的区域

### 示例

```
残基ID:   1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
共识表位: F  F  T  T  T  T  T  F  F  T  T  T  T  F  F

区域1: 3-7 (长度5) ✓
区域2: 10-13 (长度4) ✗ (< min_length=5)

最终输出: [(3, 7, avg_score, "XXXXX")]
```

## 配置参数

在 `config/default_config.yaml` 中配置：

```yaml
models:
  weights:
    bepipred: 0.25
    discotope: 0.25
    graphbepi: 0.20
    epigraph: 0.20
    vaxignml: 0.10

consensus:
  threshold: 0.5          # 共识分数阈值
  min_votes: 2            # 最小投票数
  method: "weighted_avg"  # 评分方法
```

## 输出格式

### DataFrame列说明

| 列名 | 类型 | 说明 |
|------|------|------|
| protein_id | str | 蛋白质ID |
| residue_id | int | 残基位置（1-based） |
| residue_name | str | 单字母氨基酸 |
| consensus_score | float | 共识分数 (0-1) |
| vote_count | int | 投票数 (0-5) |
| is_consensus_epitope | bool | 是否为共识表位 |
| bepipred_score | float | BepiPred分数 |
| discotope_score | float | DiscoTope分数 |
| ... | ... | 其他模型分数 |
