#!/bin/bash
# BacterialAntigenFinder 环境安装脚本
# 用于Linux系统

set -e

echo "========================================"
echo "BacterialAntigenFinder 环境安装"
echo "========================================"
echo ""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 检查conda是否安装
if ! command -v conda &> /dev/null; then
    echo -e "${RED}错误: 未找到conda，请先安装Miniconda或Anaconda${NC}"
    echo "下载地址: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo -e "${GREEN}✓ 检测到conda: $(conda --version)${NC}"
echo ""

# 获取脚本所在目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
ENVS_DIR="$PROJECT_DIR/envs"

echo "项目目录: $PROJECT_DIR"
echo "环境配置目录: $ENVS_DIR"
echo ""

# 初始化conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# 创建环境函数
create_env() {
    local env_file=$1
    local env_name=$2
    
    if [ ! -f "$env_file" ]; then
        echo -e "${YELLOW}警告: 环境文件不存在: $env_file${NC}"
        return 1
    fi
    
    echo -e "${YELLOW}创建环境: $env_name${NC}"
    
    # 检查环境是否已存在
    if conda env list | grep -q "^$env_name "; then
        echo "环境已存在，跳过创建"
        return 0
    fi
    
    # 创建环境
    if conda env create -f "$env_file" -n "$env_name"; then
        echo -e "${GREEN}✓ 环境 $env_name 创建成功${NC}"
    else
        echo -e "${RED}✗ 环境 $env_name 创建失败${NC}"
        return 1
    fi
}

# 询问是否安装所有环境
echo "可用的环境配置:"
echo "  1. master_env    - 主控环境（必需）"
echo "  2. bepipred_env  - BepiPred-3.0"
echo "  3. discotope_env - DiscoTope-3.0"
echo "  4. graphbepi_env - GraphBepi"
echo "  5. epigraph_env  - EpiGraph"
echo ""

read -p "是否安装所有环境? (y/n) " -n 1 -r
echo ""

if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo ""
    echo "开始安装所有环境..."
    echo "========================================"
    
    # 安装主控环境（必需）
    echo ""
    echo "[1/5] 安装主控环境..."
    create_env "$ENVS_DIR/master_env.yaml" "master_env"
    
    # 安装BepiPred环境
    echo ""
    echo "[2/5] 安装BepiPred-3.0环境..."
    create_env "$ENVS_DIR/bepipred_env.yaml" "bepipred_env"
    
    # 安装DiscoTope环境
    echo ""
    echo "[3/5] 安装DiscoTope-3.0环境..."
    create_env "$ENVS_DIR/discotope_env.yaml" "discotope_env"
    
    # 安装GraphBepi环境
    echo ""
    echo "[4/5] 安装GraphBepi环境..."
    create_env "$ENVS_DIR/graphbepi_env.yaml" "graphbepi_env"
    
    # 安装EpiGraph环境
    echo ""
    echo "[5/5] 安装EpiGraph环境..."
    create_env "$ENVS_DIR/epigraph_env.yaml" "epigraph_env"
    
else
    # 只安装主控环境
    echo ""
    echo "仅安装主控环境..."
    create_env "$ENVS_DIR/master_env.yaml" "master_env"
fi

echo ""
echo "========================================"
echo -e "${GREEN}环境安装完成!${NC}"
echo ""
echo "使用方法:"
echo "  1. 激活主控环境: conda activate master_env"
echo "  2. 运行Pipeline: python main.py --help"
echo ""

# 验证安装
echo "验证环境..."
conda activate master_env
python -c "import pandas; import numpy; import biopython; print('✓ 主控环境验证通过')" 2>/dev/null || echo "部分依赖可能需要手动安装"

echo ""
echo "完成!"
