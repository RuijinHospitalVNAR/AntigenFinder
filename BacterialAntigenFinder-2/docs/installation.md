# 安装指南

## 系统要求

### 硬件要求

| 组件 | 最低要求 | 推荐配置 |
|------|----------|----------|
| CPU | 4核 | 8核+ |
| 内存 | 16GB | 32GB+ |
| 存储 | 50GB | 100GB+ |
| GPU | 可选 | NVIDIA CUDA 11.3+ |

### 软件要求

| 软件 | 版本 | 说明 |
|------|------|------|
| 操作系统 | Linux (Ubuntu 18.04+) | 推荐Ubuntu 20.04 LTS |
| Python | 3.8-3.9 | 各模型要求不同 |
| Conda | Miniconda/Anaconda | 必需 |
| Docker | 19.03+ | Vaxign-ML需要 |
| CUDA | 11.3+ | GPU加速(可选) |

---

## 方式一: Conda环境安装

### 1. 安装Miniconda (如未安装)

```bash
# 下载Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 安装
bash Miniconda3-latest-Linux-x86_64.sh

# 重新加载shell
source ~/.bashrc
```

### 2. 克隆项目

```bash
git clone [repo_url] BacterialAntigenFinder
cd BacterialAntigenFinder
```

### 3. 运行安装脚本

```bash
chmod +x scripts/setup_envs.sh
./scripts/setup_envs.sh
```

安装脚本会自动创建5个Conda环境：
- `master_env` - 主控环境
- `bepipred_env` - BepiPred-3.0
- `discotope_env` - DiscoTope-3.0
- `graphbepi_env` - GraphBepi
- `epigraph_env` - EpiGraph

### 4. 手动创建环境 (可选)

如果自动安装失败，可手动创建：

```bash
# 主控环境
conda env create -f envs/master_env.yaml

# BepiPred环境
conda env create -f envs/bepipred_env.yaml

# DiscoTope环境
conda env create -f envs/discotope_env.yaml

# GraphBepi环境
conda env create -f envs/graphbepi_env.yaml

# EpiGraph环境
conda env create -f envs/epigraph_env.yaml
```

### 5. 配置模型路径

编辑 `config/default_config.yaml`:

```yaml
model_paths:
  bepipred: "/path/to/BepiPred-3.0-main"
  discotope: "/path/to/DiscoTope-3.0-master"
  graphbepi: "/path/to/GraphBepi-main"
  epigraph: "/path/to/EpiGraph-main"
  vaxignml: "/path/to/Vaxign-ML-docker-master"
```

### 6. 验证安装

```bash
conda activate master_env
python main.py --help
```

---

## 方式二: Docker安装

### 1. 安装Docker

```bash
# Ubuntu
sudo apt-get update
sudo apt-get install docker.io

# 启动Docker服务
sudo systemctl start docker
sudo systemctl enable docker

# 添加用户到docker组
sudo usermod -aG docker $USER
```

### 2. 安装NVIDIA Docker (GPU支持)

```bash
# 添加NVIDIA仓库
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
  sudo tee /etc/apt/sources.list.d/nvidia-docker.list

# 安装nvidia-docker2
sudo apt-get update
sudo apt-get install nvidia-docker2

# 重启Docker
sudo systemctl restart docker
```

### 3. 构建镜像

```bash
cd BacterialAntigenFinder

# 完整版（包含所有环境）
docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .

# 轻量版（仅主控环境）
docker build -t bacterial-antigen-finder:light -f docker/Dockerfile.light .
```

### 4. 拉取Vaxign-ML镜像

```bash
docker pull e4ong1031/vaxign-ml:latest
```

### 5. 测试运行

```bash
docker run --rm bacterial-antigen-finder:latest --help
```

---

## 模型文件准备

各预测模型需要下载模型权重文件：

### BepiPred-3.0

模型权重会在首次运行时自动从Torch Hub下载。

### DiscoTope-3.0

```bash
cd DiscoTope-3.0-master
unzip models.zip
```

### GraphBepi

模型检查点在 `model/BCE_633_GraphBepi/` 目录下。

### EpiGraph

模型检查点在 `checkpoint/` 目录下，包含10个集成模型。

### Vaxign-ML

使用Docker镜像，无需额外下载。

---

## 常见安装问题

### 问题1: Conda环境创建失败

```bash
# 更新conda
conda update conda

# 清理缓存
conda clean --all

# 重试创建
conda env create -f envs/xxx_env.yaml
```

### 问题2: PyTorch版本冲突

```bash
# 在对应环境中重新安装PyTorch
conda activate xxx_env
conda install pytorch=1.11.0 cudatoolkit=11.3 -c pytorch
```

### 问题3: Docker权限问题

```bash
# 确保用户在docker组
sudo usermod -aG docker $USER

# 重新登录生效
newgrp docker
```

### 问题4: GPU未检测到

```bash
# 验证CUDA安装
nvidia-smi

# 验证PyTorch GPU支持
python -c "import torch; print(torch.cuda.is_available())"
```

---

## 验证安装成功

### 1. 测试主控环境

```bash
conda activate master_env
python -c "from src.preprocessor import FastaParser; print('OK')"
```

### 2. 测试各预测环境

```bash
# BepiPred
conda activate bepipred_env
python -c "import torch; import esm; print('OK')"

# DiscoTope
conda activate discotope_env
python -c "import torch_geometric; print('OK')"
```

### 3. 运行测试用例

```bash
conda activate master_env
pytest tests/ -v
```

### 4. 使用示例数据测试

```bash
python main.py \
    --fasta example_data/sample_antigens.fasta \
    --pdb_dir example_data/structures/ \
    --output_dir test_output/ \
    --models bepipred
```
