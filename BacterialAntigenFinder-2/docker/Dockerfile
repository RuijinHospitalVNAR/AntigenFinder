# BacterialAntigenFinder Docker镜像
#
# 构建方式:
#   docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .
#
# 运行方式:
#   GPU_ID=0 docker compose up antigen-finder
#   # 或指定GPU:
#   GPU_ID=3 docker compose up antigen-finder

FROM continuumio/miniconda3:24.3.0-0 AS base

LABEL maintainer="BacterialAntigenFinder Team"
LABEL description="细菌抗原AI智能筛选平台"
LABEL version="1.0.0"
LABEL species.supported="Pseudomonas_aeruginosa,Klebsiella_pneumoniae,Acinetobacter_baumannii,Staphylococcus_aureus,Escherichia_coli"

# 设置环境变量
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# 安装系统依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    curl \
    dssp \
    libcifpp-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Download CCD components.cif required by mkdssp
RUN mkdir -p /var/cache/libcifpp && \
    curl --fail -L -o /var/cache/libcifpp/components.cif \
        https://files.wwpdb.org/pub/pdb/data/monomers/components.cif && \
    curl --fail -L -o /var/cache/libcifpp/mmcif_pdbx.dic \
        https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic

# 创建工作目录
WORKDIR /app

# 复制环境配置文件
COPY envs/ /app/envs/

# PyTorch CUDA 12.1 安装URL
ENV TORCH_CUDA=cu121
ENV TORCH_INDEX_URL=https://download.pytorch.org/whl/${TORCH_CUDA}

# 创建主控环境
RUN conda env create -f /app/envs/master_env.yaml && \
    conda clean -afy

# 创建BepiPred环境 + 安装CUDA 12.1版PyTorch
RUN conda env create -f /app/envs/bepipred_env.yaml && \
    /opt/conda/envs/bepipred_env/bin/pip install --no-cache-dir \
        torch==2.5.1 --index-url ${TORCH_INDEX_URL} && \
    conda clean -afy

# 创建DiscoTope环境 + 安装CUDA 12.1版PyTorch
RUN conda env create -f /app/envs/discotope_env.yaml && \
    /opt/conda/envs/discotope_env/bin/pip install --no-cache-dir \
        torch==2.5.1 --index-url ${TORCH_INDEX_URL} && \
    conda clean -afy

# 创建GraphBepi环境 + 安装CUDA 12.1版PyTorch
RUN conda env create -f /app/envs/graphbepi_env.yaml && \
    /opt/conda/envs/graphbepi_env/bin/pip install --no-cache-dir \
        torch==2.5.1 --index-url ${TORCH_INDEX_URL} && \
    conda clean -afy

# 创建EpiGraph环境 + 安装CUDA 12.1版PyTorch + torch-scatter
RUN conda env create -f /app/envs/epigraph_env.yaml && \
    /opt/conda/envs/epigraph_env/bin/pip install --no-cache-dir \
        torch==2.5.1 --index-url ${TORCH_INDEX_URL} && \
    /opt/conda/envs/epigraph_env/bin/pip install --no-cache-dir \
        torch-scatter -f https://data.pyg.org/whl/torch-2.5.1+${TORCH_CUDA}.html && \
    conda clean -afy

# 复制源代码
COPY . /app/

# 创建必要的目录（模型通过volumes挂载，不打包进镜像）
RUN mkdir -p /app/data /app/results /app/models

# 设置权限
RUN chmod +x /app/docker/entrypoint.sh 2>/dev/null || true

# 初始化conda
RUN conda init bash

# 设置默认shell
SHELL ["/bin/bash", "-c"]

# 数据和结果目录声明
VOLUME ["/app/data", "/app/results", "/app/models"]

# 暴露端口（如果需要API服务）
EXPOSE 8000

# 健康检查 - 验证conda主控环境和关键Python包可用
HEALTHCHECK --interval=60s --timeout=30s --start-period=120s --retries=3 \
    CMD conda run -n master_env python -c "import pandas; import numpy; print('Health check passed')" || exit 1

# 设置入口点 - 使用entrypoint.sh脚本
ENTRYPOINT ["/app/docker/entrypoint.sh"]

# 默认命令（显示帮助）
CMD ["--help"]
