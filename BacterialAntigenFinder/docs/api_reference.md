# API参考

本文档描述BacterialAntigenFinder的Python API，用于程序化调用。

---

## 核心类

### AntigenFinderPipeline

主Pipeline类，整合所有模块。

```python
from main import AntigenFinderPipeline

# 初始化
pipeline = AntigenFinderPipeline(config)

# 运行
results = pipeline.run(
    fasta_file="antigens.fasta",
    pdb_dir="structures/",
    organism_type="gram-",
    output_dir="results/",
    output_format="both"
)
```

#### 参数

| 参数 | 类型 | 说明 |
|------|------|------|
| config | dict | 配置字典 |

#### 方法

| 方法 | 说明 |
|------|------|
| `run(fasta_file, pdb_dir, organism_type, output_dir, output_format)` | 运行完整Pipeline |

#### 返回值

```python
{
    'candidates': pd.DataFrame,     # 候选抗原
    'consensus': pd.DataFrame,      # 共识评分
    'predictions': dict,            # 各模型预测
    'immunogenicity': pd.DataFrame  # 免疫原性评估
}
```

---

## 预处理模块

### FastaParser

```python
from src.preprocessor import FastaParser

parser = FastaParser()

# 解析FASTA文件
sequences = parser.parse("antigens.fasta")
# 返回: Dict[str, ProteinSequence]

# 获取统计信息
stats = parser.get_sequence_stats(sequences)
# 返回: {'count': 10, 'avg_length': 350, ...}

# 按长度过滤
filtered = parser.filter_by_length(sequences, min_length=50, max_length=1000)
```

### PdbValidator

```python
from src.preprocessor import PdbValidator

validator = PdbValidator()

# 验证单个PDB
structure = validator.validate("protein.pdb")
# 返回: ProteinStructure

# 验证目录
structures = validator.validate_directory("structures/")
# 返回: Dict[str, ProteinStructure]

# 提取单链
validator.extract_chain("multi_chain.pdb", "A", "chain_A.pdb")
```

### DataMapper

```python
from src.preprocessor import DataMapper

mapper = DataMapper()

# 映射序列和结构
mapped = mapper.map(sequences, structures)
# 返回: Dict[str, MappedData]

# 获取无结构的序列
no_struct = mapper.get_sequences_without_structures(mapped)
```

---

## 预测器模块

### 基类

```python
from src.predictors import BasePredictor

class BasePredictor(ABC):
    def __init__(self, model_path, env_name, threshold, use_gpu, timeout):
        ...
    
    @abstractmethod
    def predict(self, sequences, structures, organism_type) -> PredictionResult:
        ...
```

### BepipredWrapper

```python
from src.predictors import BepipredWrapper

predictor = BepipredWrapper(
    model_path="../BepiPred-3.0-main",
    env_name="bepipred_env",
    threshold=0.1512,
    use_gpu=True
)

result = predictor.predict(sequences, structures=None, organism_type="gram-")
# 返回: PredictionResult
```

### PredictionResult

```python
@dataclass
class PredictionResult:
    predictor_name: str
    predictions: List[EpitopePrediction]
    protein_scores: Dict[str, float]
    metadata: Dict[str, Any]
    
    def to_dataframe(self) -> pd.DataFrame:
        ...
    
    def get_epitopes(self, protein_id=None) -> List[EpitopePrediction]:
        ...
    
    def get_epitope_regions(self, protein_id, min_length=5) -> List[tuple]:
        ...
```

### EpitopePrediction

```python
@dataclass
class EpitopePrediction:
    protein_id: str
    residue_id: int
    residue_name: str
    score: float
    is_epitope: bool
    confidence: float = 0.0
    additional_info: Dict[str, Any] = field(default_factory=dict)
```

---

## 聚合模块

### ConsensusScorer

```python
from src.aggregator import ConsensusScorer

scorer = ConsensusScorer(
    weights={'bepipred': 0.25, 'discotope': 0.25, ...},
    threshold=0.5,
    min_votes=2,
    method='weighted_avg'
)

# 计算共识
consensus_df = scorer.compute_consensus(predictions)
# 返回: pd.DataFrame

# 获取表位区域
regions = scorer.get_epitope_regions(consensus_df, "protein_id", min_length=5)
# 返回: [(start, end, score, sequence), ...]

# 获取统计
stats = scorer.get_summary_stats(consensus_df)
```

### ImmunogenicityEvaluator

```python
from src.aggregator import ImmunogenicityEvaluator

evaluator = ImmunogenicityEvaluator(
    protegenicity_threshold=50.0,
    antigenicity_threshold=0.4
)

# 评估
immuno_df = evaluator.evaluate(vaxignml_result, sequences)
# 返回: pd.DataFrame

# 获取Top候选
top = evaluator.get_top_candidates(immuno_df, n=20, recommendation='HIGH')
```

### CandidateRanker

```python
from src.aggregator import CandidateRanker

ranker = CandidateRanker(
    top_n=50,
    min_epitope_length=5,
    min_consensus_score=0.3
)

# 排序
candidates_df = ranker.rank(consensus_df, immuno_df, sequences)
# 返回: pd.DataFrame

# 获取详情
details = ranker.get_epitope_details(candidates_df, "protein_id")
```

---

## 报告模块

### CsvExporter

```python
from src.reporter import CsvExporter

exporter = CsvExporter()

# 导出候选
path = exporter.export(candidates_df, "output/candidates.csv")

# 导出详细结果
files = exporter.export_detailed_results(consensus_df, "output/")

# 导出统计
exporter.export_summary_statistics(results, "output/stats.csv")
```

### HtmlReporter

```python
from src.reporter import HtmlReporter

reporter = HtmlReporter()

# 生成报告
path = reporter.generate(results, "output/report.html")
```

---

## 数据类

### ProteinSequence

```python
@dataclass
class ProteinSequence:
    id: str
    name: str
    description: str
    sequence: str
    length: int
    
    def to_fasta(self) -> str:
        ...
    
    def get_subsequence(self, start, end) -> str:
        ...
    
    def validate(self) -> bool:
        ...
```

### ProteinStructure

```python
@dataclass
class ProteinStructure:
    id: str
    file_path: str
    chains: List[str]
    residue_count: int
    resolution: Optional[float] = None
    is_alphafold: bool = False
    sequence: str = ""
    plddt_scores: List[float] = field(default_factory=list)
```

---

## 使用示例

### 示例1: 完整Pipeline

```python
import yaml

# 加载配置
with open('config/default_config.yaml') as f:
    config = yaml.safe_load(f)

# 修改配置
config['models']['enabled'] = ['bepipred', 'discotope']
config['consensus']['threshold'] = 0.6

# 运行Pipeline
pipeline = AntigenFinderPipeline(config)
results = pipeline.run(
    fasta_file="data/antigens.fasta",
    pdb_dir="data/structures/",
    organism_type="gram-",
    output_dir="results/"
)

# 获取候选
candidates = results['candidates']
high_priority = candidates[candidates['recommendation'] == 'HIGH']
print(f"高优先级候选数: {len(high_priority)}")
```

### 示例2: 仅使用共识评分

```python
from src.predictors import BepipredWrapper, DiscotopeWrapper
from src.aggregator import ConsensusScorer
from src.preprocessor import FastaParser, PdbValidator

# 解析数据
parser = FastaParser()
sequences = parser.parse("antigens.fasta")

validator = PdbValidator()
structures = validator.validate_directory("structures/")

# 运行预测
bp = BepipredWrapper("../BepiPred-3.0-main", "bepipred_env")
dt = DiscotopeWrapper("../DiscoTope-3.0-master", "discotope_env")

predictions = {
    'bepipred': bp.predict(sequences),
    'discotope': dt.predict(sequences, structures)
}

# 计算共识
scorer = ConsensusScorer(threshold=0.5)
consensus = scorer.compute_consensus(predictions)

# 过滤表位
epitopes = consensus[consensus['is_consensus_epitope'] == True]
print(epitopes.head(20))
```

### 示例3: 自定义权重

```python
# 增加结构基模型权重
custom_weights = {
    'bepipred': 0.15,
    'discotope': 0.35,
    'graphbepi': 0.25,
    'epigraph': 0.25,
    'vaxignml': 0.00  # 禁用
}

scorer = ConsensusScorer(weights=custom_weights, threshold=0.5)
```
