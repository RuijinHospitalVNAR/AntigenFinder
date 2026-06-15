# 预测模型详解

BacterialAntigenFinder整合了5种先进的B细胞表位预测模型。

---

## 模型概览

| 模型 | 类型 | 输入 | 方法 | 默认阈值 |
|------|------|------|------|----------|
| BepiPred-3.0 | 序列基 | FASTA | ESM-2语言模型 | 0.1512 |
| DiscoTope-3.0 | 结构基 | PDB | ESM-IF1逆折叠 | 0.90 |
| GraphBepi | 混合 | PDB/FASTA | 图神经网络 | 0.1763 |
| EpiGraph | 结构基 | PDB | 图注意力网络 | 0.1481 |
| Vaxign-ML | 序列基 | FASTA | XGBoost | 0.5 |

---

## BepiPred-3.0

### 简介

BepiPred-3.0是DTU开发的最新序列基B细胞表位预测工具，使用ESM-2蛋白质语言模型提取特征。

### 技术细节

- **输入**: 蛋白质序列 (FASTA)
- **模型**: ESM-2 (650M参数)
- **输出**: 每个残基的表位概率分数
- **阈值**: 0.1512 (可变阈值模式)

### 原理

```
序列 → ESM-2编码 → 特征提取 → 分类器 → 表位分数
```

ESM-2是Facebook AI Research开发的蛋白质语言模型，在UniRef50上预训练，能够捕获进化保守性和结构信息。

### 引用

```
Clifford JN, et al. BepiPred-3.0: Improved B-cell epitope prediction 
using protein language models. Protein Science. 2022.
```

---

## DiscoTope-3.0

### 简介

DiscoTope-3.0是DTU开发的构象表位预测工具，利用ESM-IF1逆折叠模型的潜在表示。

### 技术细节

- **输入**: 蛋白质结构 (PDB)
- **模型**: ESM-IF1 + XGBoost
- **输出**: 每个残基的表位倾向分数
- **阈值**: 0.90 (校准分数)

### 原理

```
PDB结构 → ESM-IF1嵌入 → 正向-未标记学习 → XGBoost分类 → 校准分数
```

ESM-IF1是逆折叠模型，能够从结构中推断序列，其隐藏表示包含丰富的结构-功能关联信息。

### 引用

```
Høie MH, et al. DiscoTope-3.0: improved B-cell epitope prediction 
using inverse folding latent representations. Frontiers in Immunology. 2024.
```

---

## GraphBepi

### 简介

GraphBepi是中山大学开发的图神经网络表位预测方法，结合AlphaFold预测结构和ESM-2嵌入。

### 技术细节

- **输入**: PDB结构 或 FASTA序列 (ESMFold预测)
- **模型**: Edge-enhanced GAT (EGAT)
- **输出**: 每个残基的表位分数
- **阈值**: 0.1763

### 原理

```
PDB结构 → 图构建 → ESM-2节点特征 → EGAT → 表位分类
         (残基=节点, 空间邻近=边)
```

图神经网络能够捕获残基之间的空间邻近关系，这对构象表位预测至关重要。

### 特点

- 支持ESMFold在线预测结构
- 边增强的图注意力机制

### 引用

```
Zeng Y, et al. Identifying the B-cell epitopes using AlphaFold2 
predicted structures and pretrained language model. bioRxiv. 2022.
```

---

## EpiGraph

### 简介

EpiGraph是KAIST开发的图注意力网络表位预测方法，通过捕获表位的空间聚集特性提高预测准确性。

### 技术细节

- **输入**: PDB结构
- **模型**: GAT (10折集成)
- **输出**: 每个残基的表位分数
- **阈值**: 0.1481

### 原理

```
PDB结构 → 特征提取 (ESM-2 + ESM-IF1) → 图构建 → 10个GAT模型 → 集成预测
          (RSA, 距离等)
```

10折集成提高了预测的鲁棒性和可靠性。

### 特点

- ESM-2 + ESM-IF1双模型特征
- 10折交叉验证模型集成
- 在线web服务可用

### 引用

```
Choi S, Kim D. B cell epitope prediction by capturing spatial 
clustering property of the epitopes using graph attention network. 2023.
```

---

## Vaxign-ML

### 简介

Vaxign-ML是密歇根大学开发的保护性抗原预测工具，使用机器学习方法评估蛋白质的疫苗候选潜力。

### 技术细节

- **输入**: FASTA序列 + 生物类型
- **模型**: XGBoost
- **输出**: 保护性分数 (protegenicity, 0-100百分位)
- **阈值**: 50 (百分位)

### 原理

```
序列 → 特征提取 → XGBoost分类 → 保护性分数
       ├── PSORTb (亚细胞定位)
       ├── SignalP (信号肽)
       ├── TMHMM (跨膜区)
       ├── SPAAN (粘附素)
       └── propy (理化特征)
```

### 特点

- 蛋白质级别预测（非残基级别）
- 考虑多种生物学特征
- Docker容器化部署

### 在本Pipeline中的角色

- 提供免疫原性评估
- 亚细胞定位信息
- 辅助排序（权重0.10）

### 引用

```
Ong E, et al. Vaxign-ML: Supervised Machine Learning Reverse 
Vaccinology Model for Improved Prediction of Bacterial Protective Antigens. 2020.
```

---

## 模型选择建议

### 场景1: 只有序列

```bash
--models bepipred,vaxignml
```

### 场景2: 有实验结构

```bash
--models bepipred,discotope,graphbepi,epigraph,vaxignml
```

### 场景3: 有AlphaFold结构

```bash
--models bepipred,discotope,graphbepi,epigraph,vaxignml
# DiscoTope支持AlphaFold结构，使用 --struc_type alphafold
```

### 场景4: 快速筛选

```bash
--models bepipred
```

### 场景5: 高可信度预测

```bash
--models bepipred,discotope,graphbepi,epigraph,vaxignml
--min_votes 4
--consensus_threshold 0.7
```

---

## 性能比较

基于公开基准测试的性能参考：

| 模型 | AUC | 说明 |
|------|-----|------|
| BepiPred-3.0 | ~0.65-0.70 | 序列基最优 |
| DiscoTope-3.0 | ~0.70-0.75 | 结构基最优 |
| GraphBepi | ~0.65-0.70 | GNN方法 |
| EpiGraph | ~0.65-0.70 | GAT方法 |

**注意**: 实际性能取决于数据集和评估方式。
