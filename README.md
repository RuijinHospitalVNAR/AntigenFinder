README.md
AntigenFinder
细菌抗原AI智能筛选平台 — 针对耐药细菌，整合多种B细胞表位预测模型，通过共识评分与免疫原性评估，智能输出抗原候选清单。

项目目标
耐药细菌（如铜绿假单胞菌、肺炎克雷伯菌等）的疫苗研发亟需高效的抗原筛选工具。本项目旨在构建一个具备免疫原性评估功能的AI抗原智能筛选平台，实现：

针对多种耐药细菌，自动化筛选高免疫原性抗原表位
整合多种AI预测模型，通过共识机制提升预测可靠性
输出带综合AI评分的候选抗原/表位清单，为疫苗设计提供决策支持
主要内容
支持的耐药细菌
细菌种类	中文名	耐药类型	革兰氏
Pseudomonas aeruginosa	铜绿假单胞菌	MDR/XDR	G-
Klebsiella pneumoniae	肺炎克雷伯菌	ESBL/CRE	G-
Acinetobacter baumannii	鲍曼不动杆菌	XDR	G-
Staphylococcus aureus	金黄色葡萄球菌	MRSA	G+
Escherichia coli	大肠杆菌	ESBL	G-
核心功能
多模型集成预测 — 整合5种B细胞表位预测模型
共识评分 — 加权平均 + 多数投票，降低单一模型偏差
免疫原性评估 — Vaxign-ML保护性抗原预测 + VaxiJen 2.0抗原性计算
双粒度排序 — 蛋白质级别 & 表位级别候选排序
AI综合评分 — 多维度加权评分，输出推荐等级（HIGH/MEDIUM/LOW）
Docker部署 — 完整容器化方案，支持第三方复现
实现方式
集成模型
模型	类型	输入	核心算法	权重
BepiPred-3.0	序列基	FASTA	ESM-2蛋白质语言模型	0.25
DiscoTope-3.0	结构基	PDB	ESM-IF1 + 结构特征	0.25
GraphBepi	混合	PDB	图神经网络(GNN)	0.20
EpiGraph	结构基	PDB	图注意力网络(GAT)	0.20
Vaxign-ML	序列基	FASTA	随机森林 + PSORTb	0.10
评分算法
共识评分：对每个残基，按模型权重计算加权平均分，结合多数投票（≥2个模型同意）确定表位。

免疫原性评估：

保护性抗原分数（Vaxign-ML，0-100）
抗原性分数（VaxiJen 2.0，基于氨基酸z-scale自相关系数）
亚细胞定位权重（外膜/分泌蛋白优先）
跨膜区域降权
综合评分公式：

蛋白质级别：0.35×共识分 + 0.25×保护性分 + 0.20×残基数量 + 0.20×区域数量

表位级别：0.30×共识分 + 0.25×保护性分 + 0.20×抗原性分 + 0.15×表位长度 + 0.10×模型一致性

主要流程
┌─────────────────────────────────────────────────────────────────┐
│                     AntigenFinder Pipeline                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  输入                                                            │
│  ├── FASTA 蛋白质序列文件                                        │
│  └── PDB 蛋白质结构文件                                          │
│       │                                                          │
│       ▼                                                          │
│  Step 1: 数据预处理                                              │
│  ├── FASTA序列解析 & 序列验证                                    │
│  ├── PDB结构验证                                                 │
│  └── 序列-结构映射                                               │
│       │                                                          │
│       ▼                                                          │
│  Step 2: 多模型并行预测                                          │
│  ├── BepiPred-3.0  ──→  序列基B细胞表位预测                     │
│  ├── DiscoTope-3.0  ──→  结构基构象表位预测                     │
│  ├── GraphBepi      ──→  GNN表位预测                            │
│  ├── EpiGraph       ──→  GAT表位预测                            │
│  └── Vaxign-ML      ──→  保护性抗原预测                         │
│       │                                                          │
│       ▼                                                          │
│  Step 3: 共识评分                                                │
│  ├── 加权平均 (权重: 0.25/0.25/0.20/0.20/0.10)                 │
│  ├── 多数投票 (≥2 models agree)                                 │
│  └── 共识表位识别                                                │
│       │                                                          │
│       ▼                                                          │
│  Step 4: 免疫原性评估                                            │
│  ├── Vaxign-ML 保护性抗原评分                                    │
│  ├── VaxiJen 2.0 抗原性计算 (z-scale ACC)                       │
│  ├── 亚细胞定位分析                                              │
│  └── 推荐等级判定 (HIGH/MEDIUM/LOW)                              │
│       │                                                          │
│       ▼                                                          │
│  Step 5: 候选排序                                                │
│  ├── 蛋白质级别排序 → antigen_candidates.csv                     │
│  └── 表位级别排序   → epitope_candidates.csv                     │
│       │                                                          │
│       ▼                                                          │
│  Step 6: 报告输出                                                │
│  ├── CSV 候选清单 (蛋白质级 + 表位级)                            │
│  ├── HTML 分析报告 (可视化图表)                                  │
│  ├── 详细残基级预测结果                                          │
│  └── 运行元数据 (JSON)                                           │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
安装
方式一：Conda环境安装（推荐）
# 1. 克隆项目
git clone https://github.com/RuijinHospitalVNAR/AntigenFinder.git
cd AntigenFinder/BacterialAntigenFinder

# 2. 运行环境安装脚本
chmod +x scripts/setup_envs.sh
./scripts/setup_envs.sh

# 3. 激活主控环境
conda activate master_env
方式二：Docker安装
# 轻量构建（不含模型权重，需挂载模型目录）
cd BacterialAntigenFinder
docker build -t antigen-finder:latest -f docker/Dockerfile.light .

# 完整构建（含模型，需在项目根目录执行）
docker build -t antigen-finder:latest -f docker/Dockerfile \
    --build-arg MODEL_SRC=../ .
依赖环境
Python 3.9+
PyYAML, pandas, numpy, biopython, plotly
各预测模型独立Conda环境（见 envs/ 目录）
使用方式
