# 常见问题 (FAQ)

## 安装相关

### Q1: Conda环境创建失败怎么办？

**A**: 尝试以下步骤：

```bash
# 更新conda
conda update conda

# 清理缓存
conda clean --all

# 使用更宽松的求解器
conda config --set solver libmamba
conda env create -f envs/xxx_env.yaml
```

### Q2: PyTorch CUDA版本不匹配怎么办？

**A**: 根据您的CUDA版本重新安装：

```bash
# 查看CUDA版本
nvidia-smi

# 安装对应版本
conda install pytorch=1.11.0 cudatoolkit=11.3 -c pytorch
```

### Q3: Docker镜像构建失败怎么办？

**A**: 检查网络连接，可能需要配置代理：

```bash
docker build --network=host -t bacterial-antigen-finder:latest .
```

---

## 运行相关

### Q4: 提示"模型路径不存在"怎么办？

**A**: 在 `config/default_config.yaml` 中正确配置模型路径：

```yaml
model_paths:
  bepipred: "/absolute/path/to/BepiPred-3.0-main"
  discotope: "/absolute/path/to/DiscoTope-3.0-master"
  # ...
```

### Q5: 某个模型运行失败，整个Pipeline都停止了？

**A**: 默认启用了 `--continue_on_error`，单个模型失败不会影响其他模型。确保没有禁用此选项。

### Q6: GPU内存不足怎么办？

**A**: 几种解决方案：

```bash
# 1. 使用CPU模式
python main.py --cpu_only ...

# 2. 减少每批处理的序列数
# 编辑配置文件设置batch_size

# 3. 处理较短的序列
# 过滤掉过长的序列（>1000残基）
```

### Q7: 运行速度很慢怎么办？

**A**: 优化建议：

1. 使用GPU加速
2. 只运行必要的模型: `--models bepipred,discotope`
3. 减少输入序列数量
4. 使用SSD存储

### Q8: 没有PDB结构文件怎么办？

**A**: 两种方案：

1. 只使用序列基模型：
```bash
--models bepipred,vaxignml
```

2. 使用AlphaFold预测结构：
   - 访问 https://alphafold.ebi.ac.uk/
   - 下载对应的PDB文件

---

## 结果解读

### Q9: 什么样的候选是好的候选？

**A**: 理想候选应满足：

| 指标 | 理想值 |
|------|--------|
| Recommendation | HIGH |
| Avg_Consensus_Score | ≥ 0.6 |
| Protegenicity_Score | ≥ 50 |
| Subcellular_Location | OuterMembrane/Secreted |
| Num_Epitope_Regions | ≥ 2 |

### Q10: 为什么有些蛋白质没有出现在候选列表中？

**A**: 可能原因：

1. 没有检测到符合条件的表位区域（连续≥5个残基）
2. 共识分数未达到阈值
3. 投票数不足

降低阈值可能有帮助：
```bash
--consensus_threshold 0.3 --min_votes 1
```

### Q11: 不同模型的预测结果差异很大怎么办？

**A**: 这是正常现象。不同模型使用不同的方法学：

- 序列基模型和结构基模型可能有差异
- 共识评分的设计就是为了整合这些差异
- 高一致性（vote_count ≥ 4）的预测更可靠

### Q12: 如何验证预测结果？

**A**: 建议的实验验证流程：

1. 选择Top 10-20候选
2. 合成对应的肽段
3. 进行ELISA验证
4. 用已知血清检测反应性

---

## 配置相关

### Q13: 如何调整模型权重？

**A**: 编辑 `config/default_config.yaml`：

```yaml
models:
  weights:
    bepipred: 0.30    # 增加BepiPred权重
    discotope: 0.30   # 增加DiscoTope权重
    graphbepi: 0.15
    epigraph: 0.15
    vaxignml: 0.10
```

### Q14: 如何只输出高可信度预测？

**A**: 提高阈值和投票要求：

```bash
python main.py \
    --consensus_threshold 0.7 \
    --min_votes 4 \
    ...
```

### Q15: 如何处理大量蛋白质（>100个）？

**A**: 建议分批处理：

```bash
# 拆分FASTA文件
split -l 100 large.fasta batch_

# 分批运行
for f in batch_*; do
    python main.py --fasta $f --output_dir results_$f/ ...
done

# 合并结果
cat results_*/antigen_candidates.csv > all_candidates.csv
```

---

## 数据格式

### Q16: FASTA文件格式要求？

**A**: 标准FASTA格式：

```fasta
>protein_id Description
SEQUENCE...
```

注意：
- 序列ID中避免使用特殊字符
- 序列只包含标准氨基酸字母

### Q17: PDB文件格式要求？

**A**: 
- 标准PDB格式
- 文件名应与FASTA中的序列ID匹配
- 支持多链结构

### Q18: 序列ID和PDB文件名不匹配怎么办？

**A**: 确保匹配规则：

```
FASTA: >OmpA_ECOLI ...
PDB:   OmpA_ECOLI.pdb  或  OmpA_ECOLI_A.pdb
```

程序会尝试模糊匹配。

---

## 其他问题

### Q19: 可以用于病毒蛋白吗？

**A**: 可以，设置 `--organism_type virus`。但请注意：

- Vaxign-ML针对细菌训练，病毒预测可能不够准确
- 建议主要依赖BepiPred和结构基模型

### Q20: 支持T细胞表位预测吗？

**A**: 当前版本专注于B细胞表位。T细胞表位预测需要使用其他工具如NetMHCpan。

### Q21: 如何获取更新？

**A**: 
```bash
git pull origin main
./scripts/setup_envs.sh  # 更新环境
```

### Q22: 如何报告Bug？

**A**: 在GitHub Issues中提交，包含：
- 操作系统和版本
- 错误信息
- 重现步骤
- 输入文件示例（脱敏）
