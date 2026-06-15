# 输出格式说明

## 输出文件结构

运行Pipeline后，输出目录结构如下：

```
results/
├── antigen_candidates.csv      # 候选抗原清单（主要结果）
├── analysis_report.html        # 交互式分析报告
├── run_metadata.json           # 运行元数据
└── detailed_results/           # 详细结果目录
    ├── consensus_scores.csv    # 残基级共识评分
    └── immunogenicity_scores.csv # 蛋白质级免疫原性评分
```

---

## antigen_candidates.csv

### 主要候选抗原清单

这是最重要的输出文件，包含排序后的候选抗原列表。

#### 列说明

| 列名 | 类型 | 示例 | 说明 |
|------|------|------|------|
| Rank | int | 1 | 候选排名 |
| Protein_ID | str | OmpA_ECOLI | 蛋白质ID |
| Sequence_Length | int | 346 | 序列长度 |
| Best_Epitope_Range | str | 45-60 | 最佳表位区域 |
| Best_Epitope_Sequence | str | KGFNKTGDTGVLQ | 最佳表位序列 |
| Num_Epitope_Regions | int | 3 | 表位区域数量 |
| Total_Epitope_Residues | int | 42 | 总表位残基数 |
| Avg_Consensus_Score | float | 0.7825 | 平均共识分数 |
| Max_Consensus_Score | float | 0.9234 | 最高共识分数 |
| Protegenicity_Score | float | 85.2 | 保护性分数 |
| Immunogenicity_Rank | int | 1 | 免疫原性排名 |
| Subcellular_Location | str | OuterMembrane | 亚细胞定位 |
| Recommendation | str | HIGH | 推荐等级 |
| Composite_Score | float | 0.8123 | 综合评分 |
| All_Epitope_Regions | str | 45-60(0.78); 120-135(0.65) | 所有表位区域 |

#### 示例内容

```csv
Rank,Protein_ID,Sequence_Length,Best_Epitope_Range,Best_Epitope_Sequence,Num_Epitope_Regions,Total_Epitope_Residues,Avg_Consensus_Score,Max_Consensus_Score,Protegenicity_Score,Immunogenicity_Rank,Subcellular_Location,Recommendation,Composite_Score,All_Epitope_Regions
1,OmpA_ECOLI,346,45-60,KGFNKTGDTGVLQ,3,42,0.7825,0.9234,85.2,1,OuterMembrane,HIGH,0.8123,45-60(0.78); 120-135(0.65); 200-215(0.72)
2,FimH_ECOLI,300,120-138,NNPVTGQGTANVYV,2,28,0.7112,0.8456,78.5,3,Fimbrial,HIGH,0.7456,120-138(0.71); 245-258(0.68)
```

---

## consensus_scores.csv

### 残基级共识评分结果

包含每个蛋白质每个残基的预测详情。

#### 列说明

| 列名 | 类型 | 说明 |
|------|------|------|
| protein_id | str | 蛋白质ID |
| residue_id | int | 残基位置（1-based） |
| residue_name | str | 单字母氨基酸 |
| consensus_score | float | 共识分数 (0-1) |
| vote_count | int | 投票数 (0-5) |
| is_consensus_epitope | bool | 是否为共识表位 |
| bepipred_score | float | BepiPred分数 |
| bepipred_epitope | bool | BepiPred预测结果 |
| discotope_score | float | DiscoTope分数 |
| discotope_epitope | bool | DiscoTope预测结果 |
| graphbepi_score | float | GraphBepi分数 |
| graphbepi_epitope | bool | GraphBepi预测结果 |
| epigraph_score | float | EpiGraph分数 |
| epigraph_epitope | bool | EpiGraph预测结果 |
| vaxignml_score | float | Vaxign-ML分数 |
| vaxignml_epitope | bool | Vaxign-ML预测结果 |

#### 示例内容

```csv
protein_id,residue_id,residue_name,consensus_score,vote_count,is_consensus_epitope,bepipred_score,bepipred_epitope,discotope_score,discotope_epitope,...
OmpA_ECOLI,1,M,0.25,1,False,0.12,False,0.35,False,...
OmpA_ECOLI,2,K,0.32,1,False,0.28,False,0.40,False,...
OmpA_ECOLI,45,K,0.82,4,True,0.85,True,0.78,True,...
```

---

## immunogenicity_scores.csv

### 蛋白质级免疫原性评分

包含每个蛋白质的免疫原性评估结果。

#### 列说明

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
| recommendation | str | 推荐等级 (HIGH/MEDIUM/LOW) |
| composite_score | float | 综合评分 |

---

## analysis_report.html

### 交互式分析报告

HTML格式的可视化报告，包含：

1. **摘要统计卡片**
   - 分析蛋白质数
   - 总残基数
   - 表位残基数
   - 高优先级候选数

2. **候选抗原表格**
   - Top 20候选详情
   - 可排序列
   - 推荐等级颜色标识

3. **可视化图表**
   - 推荐分布饼图
   - 共识分数直方图

4. **预测器统计**
   - 各模型预测表位数量
   - 模型间一致性

---

## run_metadata.json

### 运行元数据

记录Pipeline运行的详细信息。

```json
{
  "version": "1.0.0",
  "run_time": "2024-01-15T10:30:45.123456",
  "elapsed_seconds": 325.67,
  "input": {
    "fasta_file": "/path/to/antigens.fasta",
    "pdb_dir": "/path/to/structures/",
    "organism_type": "gram-",
    "num_sequences": 50,
    "num_structures": 45
  },
  "models": [
    "bepipred",
    "discotope",
    "graphbepi",
    "epigraph",
    "vaxignml"
  ],
  "output": {
    "num_candidates": 35,
    "output_dir": "/path/to/results/"
  }
}
```

---

## 解读指南

### 如何选择候选抗原

1. **查看Recommendation列**
   - `HIGH`: 优先考虑
   - `MEDIUM`: 次选
   - `LOW`: 不推荐

2. **评估表位质量**
   - `Avg_Consensus_Score` ≥ 0.6 为优质表位
   - `Num_Epitope_Regions` ≥ 2 表示多个潜在靶点

3. **考虑亚细胞定位**
   - 优选: OuterMembrane, Secreted, Fimbrial
   - 避免: Cytoplasmic

4. **综合评分排序**
   - `Composite_Score` 是综合考虑所有因素的最终评分
   - 建议选择Top 10-20进行实验验证
