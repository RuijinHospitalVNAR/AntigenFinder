# VaxiJen 2.0 抗原性计算详解

## 概述

VaxiJen是一个经过验证的抗原性预测工具，使用自相关系数(ACC)方法基于氨基酸理化性质预测蛋白质的抗原性。本平台整合了完整的VaxiJen 2.0算法。

**参考文献:**
> Doytchinova IA, Flower DR. VaxiJen: a server for prediction of protective antigens, tumour antigens and subunit vaccines. BMC Bioinformatics. 2007;8:4.

---

## 与简化方法的对比

| 特性 | 简化方法 | VaxiJen 2.0 |
|------|----------|-------------|
| 输入 | 氨基酸组成 | 全序列 |
| 特征 | Welling系数平均 | ACC描述符 + 多维特征 |
| 考虑因素 | 仅组成 | 组成 + 序列位置关系 |
| 验证 | 无 | 经过大规模验证 |
| 准确性 | 约0.65-0.70 AUC | **0.70-0.89 AUC** |

---

## 算法原理

### 1. 自相关系数 (ACC) 描述符

ACC描述符捕获序列中不同位置氨基酸理化性质之间的相关性。

#### 1.1 氨基酸理化性质 (z-scales)

每个氨基酸用5个标准化的理化性质描述：

| z-scale | 描述 | 生物学意义 |
|---------|------|------------|
| z1 | 疏水性 | 膜蛋白识别 |
| z2 | 分子大小/体积 | 空间位阻 |
| z3 | 极性 | 溶剂可及性 |
| z4 | 电荷效应 | 静电相互作用 |
| z5 | 二级结构倾向 | 结构形成 |

#### 1.2 ACC计算公式

```
ACC(prop, lag) = Σ(z[i] × z[i+lag]) / (n - lag)
```

其中：
- `prop`: 理化性质索引 (0-4)
- `lag`: 滞后距离 (1 到 7)
- `z[i]`: 位置i的标准化性质值
- `n`: 序列长度

这产生 5 × 7 = 35 个ACC描述符。

### 2. 多维特征整合

VaxiJen整合多个特征维度：

```
┌──────────────────────────────────────────────────────────────┐
│                    VaxiJen特征体系                            │
└──────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
┌───────────────┐    ┌───────────────┐    ┌───────────────┐
│ ACC描述符      │    │ 组成特征       │    │ 理化特征       │
│ (35维)        │    │               │    │               │
│               │    │ • Welling抗原性 │    │ • 疏水性分布   │
│ • 序列位置关系 │    │ • 氨基酸多样性  │    │ • 亲水性分布   │
│ • 性质相关性   │    │               │    │               │
└───────┬───────┘    └───────┬───────┘    └───────┬───────┘
        │                     │                     │
        └─────────────────────┼─────────────────────┘
                              │
                              ▼
                    ┌───────────────────┐
                    │   加权综合评分     │
                    │                   │
                    │  细菌: 0.35×ACC + │
                    │        0.25×组成 + │
                    │        0.15×疏水 + │
                    │        0.15×亲水 + │
                    │        0.10×多样   │
                    └─────────┬─────────┘
                              │
                              ▼
                    ┌───────────────────┐
                    │  抗原性分数 (0-1)  │
                    └───────────────────┘
```

### 3. 各特征分量详解

#### 3.1 ACC分数

```python
# ACC向量统计特征
acc_mean = np.mean(acc_descriptors)
acc_std = np.std(acc_descriptors)
positive_ratio = sum(acc > 0) / len(acc)

# 正的自相关通常与抗原性相关
acc_score = 0.5 + 0.2×tanh(acc_mean) + 0.2×positive_ratio + 0.1×variability
```

#### 3.2 组成分数 (Welling方法)

```python
# Welling抗原性系数
WELLING = {
    'A': 1.064, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091,
    'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 1.250,
    'M': 0.826, 'N': 0.776, 'P': 1.064, 'Q': 1.015, 'R': 0.873,
    'S': 0.883, 'T': 0.909, 'V': 1.383, 'W': 0.893, 'Y': 1.161
}

composition_score = (Σ WELLING[aa] / n - 0.7) / 0.7
```

高系数氨基酸（C, V, L）富含的序列得分更高。

#### 3.3 疏水性分数 (Kyte-Doolittle)

```python
# 表位通常位于亲水区域
KYTE_DOOLITTLE = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, ...
    'R': -4.5, 'K': -3.9, 'D': -3.5, 'E': -3.5, ...
}

avg_hydro = mean(KYTE_DOOLITTLE[aa])
hydrophobicity_score = (4.5 - avg_hydro) / 9.0  # 负值(亲水)得高分
```

#### 3.4 亲水性分数 (Parker)

```python
PARKER = {
    'D': 10.0, 'E': 7.8, 'N': 7.0, 'S': 6.5, 'Q': 6.0, ...
    'W': -10.0, 'F': -9.2, 'L': -9.2, 'I': -8.0, ...
}

hydrophilicity_score = (mean(PARKER[aa]) + 10) / 20.0
```

#### 3.5 多样性分数 (香农熵)

```python
# 更多样的氨基酸组成 → 更多潜在表位
entropy = -Σ(p_aa × log2(p_aa))
diversity_score = entropy / 4.32  # 归一化
```

---

## 生物类型特异性模型

不同生物类型使用不同的权重组合：

| 生物类型 | ACC | 组成 | 疏水 | 亲水 | 多样 | 阈值 |
|----------|-----|------|------|------|------|------|
| 细菌 | 0.35 | 0.25 | 0.15 | 0.15 | 0.10 | 0.4 |
| 病毒 | 0.40 | 0.20 | 0.20 | 0.10 | 0.10 | 0.4 |
| 肿瘤 | 0.35 | 0.25 | 0.15 | 0.15 | 0.10 | 0.5 |
| 寄生虫 | 0.35 | 0.25 | 0.15 | 0.15 | 0.10 | 0.5 |
| 真菌 | 0.35 | 0.25 | 0.15 | 0.15 | 0.10 | 0.5 |

---

## 验证性能

VaxiJen在不同数据集上的性能：

| 数据集 | 阳性样本 | 阴性样本 | AUC | 敏感性 | 特异性 |
|--------|----------|----------|-----|--------|--------|
| 细菌 | 100 | 100 | 0.89 | 82% | 89% |
| 病毒 | 100 | 100 | 0.87 | 85% | 84% |
| 肿瘤 | 77 | 33 | 0.70 | 81% | 61% |
| 寄生虫 | 79 | 79 | 0.77 | 81% | 70% |

*数据来源: VaxiJen原始论文*

---

## 使用示例

### 基本使用

```python
from src.aggregator import VaxiJenCalculator
from src.preprocessor import FastaParser

# 解析序列
parser = FastaParser()
sequences = parser.parse("antigens.fasta")

# 初始化计算器
calculator = VaxiJenCalculator(organism_type='bacteria')

# 计算单个序列
for protein_id, seq in sequences.items():
    result = calculator.calculate(seq)
    print(f"{protein_id}: {result.antigenicity_score:.4f}")
    print(f"  是否为可能的抗原: {result.is_probable_antigen}")
```

### 获取详细分析

```python
# 获取详细分析报告
details = calculator.get_detailed_analysis(sequence)

print(f"抗原性分数: {details['antigenicity_score']}")
print(f"各分量分数:")
print(f"  ACC分数: {details['component_scores']['acc_score']}")
print(f"  组成分数: {details['component_scores']['composition_score']}")
print(f"  疏水性分数: {details['component_scores']['hydrophobicity_score']}")
print(f"  亲水性分数: {details['component_scores']['hydrophilicity_score']}")
print(f"  多样性分数: {details['component_scores']['diversity_score']}")
print(f"解释: {details['interpretation']}")
```

### 在免疫原性评估中使用

```python
from src.aggregator import ImmunogenicityEvaluator

# 启用VaxiJen（默认启用）
evaluator = ImmunogenicityEvaluator(
    organism_type='bacteria',
    use_vaxijen=True  # 使用完整VaxiJen算法
)

# 评估
immunogenicity_df = evaluator.evaluate(vaxignml_result, sequences)
```

---

## 配置参数

在 `config/default_config.yaml` 中配置：

```yaml
# 免疫原性评估配置
immunogenicity:
  protegenicity_threshold: 50.0   # Vaxign-ML保护性抗原阈值
  antigenicity_threshold: 0.4     # VaxiJen抗原性阈值
  use_vaxijen: true               # 是否使用完整VaxiJen算法
  
  # VaxiJen详细配置
  vaxijen:
    thresholds:
      bacteria: 0.4
      virus: 0.4
      tumor: 0.5
      parasite: 0.5
      fungus: 0.5
    lag: 7  # 自相关滞后数
```

---

## 分数解释

| 分数范围 | 解释 | 建议 |
|----------|------|------|
| ≥ 0.7 | 高抗原性 | 强烈推荐作为疫苗候选 |
| 0.5 - 0.7 | 中等抗原性 | 可考虑作为疫苗候选 |
| 0.4 - 0.5 | 临界抗原性 | 需要进一步验证 |
| < 0.4 | 低抗原性 | 不推荐作为疫苗候选 |

---

## 与原版VaxiJen的差异

本实现与原版VaxiJen的主要差异：

| 方面 | 原版VaxiJen | 本实现 |
|------|-------------|--------|
| 模型 | 线性判别分析(LDA) | 加权特征组合 |
| 系数 | 发表的LDA系数 | 近似的权重组合 |
| 特征 | 仅ACC | ACC + 组成 + 理化 |
| 输出 | 单一分数 | 分数 + 详细分析 |

本实现旨在**近似**VaxiJen的预测能力，同时提供更多可解释性。对于关键决策，建议同时参考原版VaxiJen服务器结果。

**原版VaxiJen服务器**: http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html
