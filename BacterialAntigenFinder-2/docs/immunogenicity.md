# 免疫原性评估详解

## 概述

免疫原性评估模块(`ImmunogenicityEvaluator`)整合Vaxign-ML预测结果，结合亚细胞定位、抗原性评分等多个维度，对蛋白质进行免疫原性排名。

## 评估维度

### 1. 保护性分数 (Protegenicity Score)

**来源**: Vaxign-ML预测

**含义**: 该蛋白质作为保护性抗原的可能性，以百分位数表示(0-100)。

**阈值**: 默认 ≥50 判定为保护性抗原

```python
is_protective_antigen = protegenicity_score >= 50.0
```

### 2. 亚细胞定位 (Subcellular Location)

**来源**: Vaxign-ML内置的PSORTb预测

**优选位置**: 表面暴露的蛋白质更适合作为疫苗靶标

| 定位 | 权重 | 说明 |
|------|------|------|
| OuterMembrane | 1.5 | 外膜蛋白，直接暴露于免疫系统 |
| Fimbrial | 1.5 | 菌毛蛋白，表面暴露 |
| Surface | 1.5 | 表面蛋白 |
| Extracellular | 1.4 | 胞外/分泌蛋白 |
| Secreted | 1.4 | 分泌蛋白 |
| Periplasmic | 1.3 | 周质蛋白（G-菌） |
| CellWall | 1.3 | 细胞壁蛋白（G+菌） |
| InnerMembrane | 1.0 | 内膜蛋白（中性） |
| Unknown | 1.0 | 未知（中性） |
| Cytoplasmic | 0.8 | 胞质蛋白（降权） |

### 3. 抗原性分数 (Antigenicity Score)

**来源**: 基于氨基酸组成计算

**方法**: 类似VaxiJen的自相关特征方法（简化版）

每个氨基酸有一个抗原性贡献系数：

```python
antigenicity_coefs = {
    'A': 1.064, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091,
    'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 1.250,
    'M': 0.826, 'N': 0.776, 'P': 1.064, 'Q': 1.015, 'R': 0.873,
    'S': 0.883, 'T': 0.909, 'V': 1.383, 'W': 0.893, 'Y': 1.161
}
```

**计算**:
```python
antigenicity_score = Σ(coef_aa) / sequence_length
```

系数 > 1 表示该氨基酸有利于抗原性（如 Cys=1.412, Val=1.383）。

### 4. 跨膜区域 (Transmembrane Regions)

**来源**: Vaxign-ML内置的TMHMM预测

**惩罚规则**: 过多跨膜区域的蛋白难以表达和纯化

```python
tm_penalty = 0.8 if transmembrane_regions > 2 else 1.0
```

## 综合评分公式

```python
composite_score = protegenicity_score × location_weight × antigenicity_score × tm_penalty
```

### 示例计算

| 指标 | 值 | 说明 |
|------|-----|------|
| protegenicity_score | 75.0 | Vaxign-ML预测 |
| subcellular_location | OuterMembrane | → weight=1.5 |
| antigenicity_score | 1.05 | 序列组成计算 |
| transmembrane_regions | 1 | → penalty=1.0 |

```
composite_score = 75.0 × 1.5 × 1.05 × 1.0 = 118.125
```

## 推荐等级判定

### 判定逻辑

```python
if is_protective_antigen AND location in PREFERRED_LOCATIONS AND antigenicity >= 0.4:
    recommendation = 'HIGH'
elif is_protective_antigen OR location in PREFERRED_LOCATIONS:
    recommendation = 'MEDIUM'
else:
    recommendation = 'LOW'
```

### 等级定义

| 等级 | 条件 | 建议 |
|------|------|------|
| **HIGH** | 保护性抗原 + 优选定位 + 高抗原性 | 优先进行实验验证 |
| **MEDIUM** | 满足保护性或优选定位之一 | 次优候选，可考虑 |
| **LOW** | 均不满足 | 不推荐作为疫苗靶标 |

### 优选亚细胞定位

```python
PREFERRED_LOCATIONS = [
    'OuterMembrane',
    'Extracellular',
    'Periplasmic',
    'CellWall',
    'Secreted',
    'Fimbrial',
    'Surface'
]
```

## 评估流程

```
┌─────────────────────────────────────────────────────────────┐
│  输入: Vaxign-ML预测结果 + 蛋白质序列                         │
└─────────────────────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 1: 提取蛋白质级别信息                                   │
│  • protegenicity_score (保护性分数)                          │
│  • subcellular_location (亚细胞定位)                         │
│  • has_signal_peptide (信号肽)                               │
│  • transmembrane_regions (跨膜区数量)                        │
└─────────────────────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 2: 计算抗原性分数                                       │
│  • 统计序列中各氨基酸组成                                     │
│  • 计算抗原性系数加权平均                                     │
└─────────────────────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 3: 计算综合评分                                         │
│  composite = protegenicity × location_weight ×               │
│              antigenicity × tm_penalty                       │
└─────────────────────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 4: 排名和推荐                                           │
│  • 按综合评分排序                                             │
│  • 分配排名 (immunogenicity_rank)                            │
│  • 判定推荐等级 (HIGH/MEDIUM/LOW)                            │
└─────────────────────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────┐
│  输出: ImmunogenicityResult DataFrame                        │
└─────────────────────────────────────────────────────────────┘
```

## 输出格式

### DataFrame列说明

| 列名 | 类型 | 说明 |
|------|------|------|
| protein_id | str | 蛋白质ID |
| protegenicity_score | float | Vaxign-ML保护性分数 (0-100) |
| is_protective_antigen | bool | 是否为保护性抗原 |
| subcellular_location | str | 亚细胞定位 |
| has_signal_peptide | bool | 是否有信号肽 |
| transmembrane_regions | int | 跨膜区域数量 |
| antigenicity_score | float | 抗原性分数 |
| immunogenicity_rank | int | 免疫原性排名 |
| recommendation | str | 推荐等级 |
| composite_score | float | 综合评分 |

## 与共识评分的整合

最终候选排序综合考虑：

1. **表位质量**: 共识评分 (残基级别)
2. **蛋白质质量**: 免疫原性评估 (蛋白质级别)

```python
final_score = 0.35 × consensus_component +
              0.25 × protegenicity_component +
              0.20 × epitope_residue_component +
              0.20 × epitope_region_component
```

## 配置参数

```yaml
# 在配置文件中可调整
immunogenicity:
  protegenicity_threshold: 50.0   # 保护性抗原阈值
  antigenicity_threshold: 0.4     # 抗原性阈值
```
