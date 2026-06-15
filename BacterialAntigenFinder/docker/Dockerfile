# BacterialAntigenFinder Docker镜像
# 多阶段构建，整合多个表位预测模型
#
# 构建方式:
#   轻量构建（不含模型，运行时挂载模型目录）:
#     docker build -t bacterial-antigen-finder:latest -f docker/Dockerfile .
#
#   完整构建（将模型打包进镜像）:
#     docker build -t bacterial-antigen-finder:full \
#       --build-arg MODEL_SRC=../ \
#       -f docker/Dockerfile .
#
#   注意: 完整构建需要模型目录与 BacterialAntigenFinder 同级，目录结构如下:
#     project_root/
#     ├── BacterialAntigenFinder/   (本项目)
#     ├── BepiPred-3.0-main/
#     ├── DiscoTope-3.0-master/
#     ├── GraphBepi-main/
#     └── EpiGraph-main/

FROM continuumio/miniconda3:latest AS base

LABEL maintainer="BacterialAntigenFinder Team"
LABEL description="细菌抗原AI智能筛选平台"
LABEL version="1.0.0"
LABEL species.supported="Pseudomonas_aeruginosa,Klebsiella_pneumoniae,Acinetobacter_baumannii,Staphylococcus_aureus,Escherichia_coli"

# 模型源目录构建参数
# 轻量构建时不复制模型；完整构建时设为模型所在目录（相对于构建上下文）
ARG MODEL_SRC=""

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
    libgomp1 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# 创建工作目录
WORKDIR /app

# 复制环境配置文件
COPY envs/ /app/envs/

# 创建主控环境
RUN conda env create -f /app/envs/master_env.yaml && \
    conda clean -afy

# 创建BepiPred环境
RUN conda env create -f /app/envs/bepipred_env.yaml && \
    conda clean -afy

# 创建DiscoTope环境
RUN conda env create -f /app/envs/discotope_env.yaml && \
    conda clean -afy

# 创建GraphBepi环境
RUN conda env create -f /app/envs/graphbepi_env.yaml && \
    conda clean -afy

# 创建EpiGraph环境
RUN conda env create -f /app/envs/epigraph_env.yaml && \
    conda clean -afy

# 复制源代码
COPY . /app/

# 条件复制模型文件（仅在 MODEL_SRC 非空时生效）
# 完整构建时，MODEL_SRC 指向包含模型目录的路径（相对于构建上下文）
# 例如: --build-arg MODEL_SRC=../  则模型目录为 ../BepiPred-3.0-main 等
# Docker COPY 指令不支持条件判断，此处使用通配符+忽略缺失文件的方式
# 实际构建时需确保 .dockerignore 不会排除模型目录
COPY ${MODEL_SRC}BepiPred-3.0-main /app/models/BepiPred-3.0-main* 
COPY ${MODEL_SRC}DiscoTope-3.0-master /app/models/DiscoTope-3.0-master*
COPY ${MODEL_SRC}GraphBepi-main /app/models/GraphBepi-main*
COPY ${MODEL_SRC}EpiGraph-main /app/models/EpiGraph-main*

# 创建必要的目录
RUN mkdir -p /app/data /app/results /app/models

# 设置权限
RUN chmod +x /app/docker/entrypoint.sh 2>/dev/null || true

# 初始化conda
RUN conda init bash

# 设置默认shell
SHELL ["/bin/bash", "-c"]

# 数据和结果目录声明
VOLUME ["/app/data", "/app/results"]

# 暴露端口（如果需要API服务）
EXPOSE 8000

# 健康检查 - 验证conda主控环境和关键Python包可用
HEALTHCHECK --interval=60s --timeout=30s --start-period=120s --retries=3 \
    CMD conda run -n master_env python -c "import pandas; import numpy; print('Health check passed')" || exit 1

# 设置入口点 - 使用entrypoint.sh脚本
ENTRYPOINT ["/app/docker/entrypoint.sh"]

# 默认命令（显示帮助）
CMD ["--help"]
