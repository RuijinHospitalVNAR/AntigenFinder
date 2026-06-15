# BacterialAntigenFinder 文档

## 细菌抗原AI智能筛选平台

**版本**: 1.0.0  
**最后更新**: 2024年

---

## 目录

1. [概述](#概述)
2. [系统架构](architecture.md)
3. [安装指南](installation.md)
4. [快速入门](quickstart.md)
5. [命令行参数](cli_reference.md)
6. [模型详解](models.md)
7. [共识评分算法](consensus_algorithm.md)
8. [免疫原性评估](immunogenicity.md)
9. [**VaxiJen抗原性计算**](vaxijen_algorithm.md) ⭐ 新增
10. [**综合推荐流程**](recommendation_pipeline.md) ⭐
11. [输出格式说明](output_format.md)
11. [API参考](api_reference.md)
12. [常见问题](faq.md)

---

## 概述

**BacterialAntigenFinder** 是一个专门针对耐药细菌的抗原智能筛选平台，通过整合多种B细胞表位预测模型，结合共识评分和免疫原性评估，输出高质量的疫苗候选抗原清单。

### 核心功能

| 功能模块 | 描述 |
|----------|------|
| **多模型集成** | 整合5种先进的表位预测模型 |
| **共识评分** | 加权平均 + 多数投票双重验证 |
| **免疫原性评估** | 基于Vaxign-ML的保护性抗原预测 |
| **智能排序** | 多维度综合评分排序 |
| **灵活部署** | 支持Conda和Docker两种部署方式 |

### 支持的预测模型

| 模型 | 类型 | 方法 | 输入 |
|------|------|------|------|
| BepiPred-3.0 | 序列基 | ESM-2蛋白质语言模型 | FASTA |
| DiscoTope-3.0 | 结构基 | ESM-IF1逆折叠模型 | PDB |
| GraphBepi | 混合 | 图神经网络(GNN) | PDB/FASTA |
| EpiGraph | 结构基 | 图注意力网络(GAT) | PDB |
| Vaxign-ML | 序列基 | XGBoost机器学习 | FASTA |

### 适用场景

- 细菌疫苗候选抗原筛选
- 耐药菌株保护性抗原发现
- B细胞表位预测和验证
- 免疫原性快速评估

---

## 快速链接

- [GitHub仓库](#)
- [问题反馈](#)
- [更新日志](changelog.md)
