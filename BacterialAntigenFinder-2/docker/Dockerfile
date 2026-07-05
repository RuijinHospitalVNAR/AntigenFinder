# BacterialAntigenFinder Docker镜像
#
# 构建方式:
#   docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .
#
# 运行方式:
#   GPU_ID=0 docker compose up antigen-finder
#   # 或指定GPU:
#   GPU_ID=3 docker compose up antigen-finder

FROM docker.m.daocloud.io/continuumio/miniconda3:latest AS base

LABEL maintainer="BacterialAntigenFinder Team"
LABEL description="细菌抗原AI智能筛选平台"
LABEL version="1.0.0"
LABEL species.supported="Pseudomonas_aeruginosa,Klebsiella_pneumoniae,Acinetobacter_baumannii,Staphylococcus_aureus,Escherichia_coli"

# 设置环境变量
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
# 禁用conda repodata分片下载（网络不稳定时容易失败）
ENV CONDA_REPODATA_USE_SHARDING=false

# 替换为国内镜像源（国内网络加速）
RUN conda config --set show_channel_urls yes && \
    conda config --set remote_max_retries 10 && \
    conda config --set remote_read_timeout_secs 120 && \
    conda config --remove-key channels 2>/dev/null; \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ && \
    printf 'default_channels:\n  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main\n  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free\n  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r\noverride_channels_enabled: true\n' >> /opt/conda/.condarc && \
    sed -i 's|deb.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list.d/debian.sources

# 设置pip国内镜像 + 重试
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
ENV PIP_RETRIES=10
ENV PIP_TIMEOUT=120

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
# (Skipped in build: only needed for PDB structure analysis.
#  To enable, uncomment below or download at runtime.)
RUN mkdir -p /var/cache/libcifpp && \
    touch /var/cache/libcifpp/components.cif /var/cache/libcifpp/mmcif_pdbx.dic

# 创建工作目录
WORKDIR /app

# 复制环境配置文件
COPY envs/ /app/envs/

# PyTorch CUDA 12.1 wheels (交大镜像)
ARG PYTORCH_MIRROR=https://mirror.sjtu.edu.cn/pytorch-wheels/cu121

# Helper: pip install with retry for unstable connections
RUN echo '#!/bin/bash\nretry_pip() {\n  local env=$1; shift\n  for i in 1 2 3 4 5 6 7 8 9 10; do\n    echo "[retry_pip] Attempt $i/10: /opt/conda/envs/$env/bin/pip install $@" >&2\n    if /opt/conda/envs/$env/bin/pip install --no-cache-dir "$@"; then\n      return 0\n    fi\n    echo "[retry_pip] Failed, waiting 10s..." >&2\n    sleep 10\n  done\n  return 1\n}\nretry_pip "$@"' > /usr/local/bin/retry_pip && chmod +x /usr/local/bin/retry_pip

# Helper: conda env create with retry (HTTP errors are intermittent)
RUN echo '#!/bin/bash\nfor i in 1 2 3 4 5 6 7 8 9 10; do\n  echo "[retry_conda] Attempt $i/10: conda env create -f $@" >&2\n  if conda env create -f "$@"; then\n    conda clean -afy\n    exit 0\n  fi\n  echo "[retry_conda] Failed, removing partial env and waiting 15s..." >&2\n  sleep 15\n  # 获取环境名并删除可能的部分环境\n  env_name=$(grep -m1 "^name:" "$@" | awk "{print \$2}")\n  [ -n "$env_name" ] && conda env remove -n "$env_name" -y 2>/dev/null\n  rm -rf /opt/conda/envs/$env_name 2>/dev/null\ndone\nexit 1' > /usr/local/bin/retry_conda && chmod +x /usr/local/bin/retry_conda

# 创建主控环境
RUN retry_conda /app/envs/master_env.yaml

# 创建BepiPred环境
RUN retry_conda /app/envs/bepipred_env.yaml

# BepiPred: 安装PyTorch (重试10次)
RUN retry_pip bepipred_env torch==2.5.1 -f ${PYTORCH_MIRROR}

# 创建DiscoTope环境
RUN retry_conda /app/envs/discotope_env.yaml

# DiscoTope: 安装PyTorch
RUN retry_pip discotope_env torch==2.5.1 -f ${PYTORCH_MIRROR}

# 创建GraphBepi环境
RUN retry_conda /app/envs/graphbepi_env.yaml

# GraphBepi: 安装PyTorch
RUN retry_pip graphbepi_env torch==2.5.1 -f ${PYTORCH_MIRROR}

# 创建EpiGraph环境
RUN retry_conda /app/envs/epigraph_env.yaml

# EpiGraph: 安装PyTorch + torch-scatter
RUN retry_pip epigraph_env torch==2.5.1 -f ${PYTORCH_MIRROR} && \
    /opt/conda/envs/epigraph_env/bin/pip install --no-cache-dir \
        torch-scatter -f https://data.pyg.org/whl/torch-2.5.1+cu121.html

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
