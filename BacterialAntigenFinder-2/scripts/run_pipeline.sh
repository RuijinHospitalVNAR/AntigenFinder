#!/bin/bash
# BacterialAntigenFinder 运行脚本
# 用于Linux系统

set -e

echo "========================================"
echo "BacterialAntigenFinder - 抗原筛选Pipeline"
echo "========================================"
echo ""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# 获取脚本目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# 默认参数
FASTA_FILE=""
PDB_DIR=""
OUTPUT_DIR="./results"
ORGANISM_TYPE="gram-"
MODELS="bepipred,discotope,graphbepi,epigraph,vaxignml"
USE_GPU=true
LOG_LEVEL="INFO"

# 显示帮助
show_help() {
    echo "用法: $0 [选项]"
    echo ""
    echo "必需参数:"
    echo "  -f, --fasta FILE       FASTA序列文件"
    echo "  -p, --pdb_dir DIR      PDB结构目录"
    echo ""
    echo "可选参数:"
    echo "  -o, --output DIR       输出目录 (默认: ./results)"
    echo "  -t, --organism TYPE    生物类型: gram+, gram-, virus (默认: gram-)"
    echo "  -m, --models MODELS    启用的模型，逗号分隔 (默认: 全部)"
    echo "  --cpu                  仅使用CPU"
    echo "  --log LEVEL            日志级别: DEBUG, INFO, WARNING, ERROR"
    echo "  -h, --help             显示帮助"
    echo ""
    echo "示例:"
    echo "  $0 -f antigens.fasta -p structures/ -o results/"
    echo "  $0 -f antigens.fasta -p structures/ -t gram+ --cpu"
    echo ""
}

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--fasta)
            FASTA_FILE="$2"
            shift 2
            ;;
        -p|--pdb_dir)
            PDB_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--organism)
            ORGANISM_TYPE="$2"
            shift 2
            ;;
        -m|--models)
            MODELS="$2"
            shift 2
            ;;
        --cpu)
            USE_GPU=false
            shift
            ;;
        --log)
            LOG_LEVEL="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo -e "${RED}未知参数: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# 验证必需参数
if [ -z "$FASTA_FILE" ]; then
    echo -e "${RED}错误: 必须指定FASTA文件 (-f)${NC}"
    show_help
    exit 1
fi

if [ -z "$PDB_DIR" ]; then
    echo -e "${RED}错误: 必须指定PDB目录 (-p)${NC}"
    show_help
    exit 1
fi

# 验证文件存在
if [ ! -f "$FASTA_FILE" ]; then
    echo -e "${RED}错误: FASTA文件不存在: $FASTA_FILE${NC}"
    exit 1
fi

if [ ! -d "$PDB_DIR" ]; then
    echo -e "${RED}错误: PDB目录不存在: $PDB_DIR${NC}"
    exit 1
fi

# 初始化conda
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true

# 激活环境
echo "激活主控环境..."
conda activate master_env 2>/dev/null || {
    echo -e "${YELLOW}警告: 无法激活master_env，尝试直接运行${NC}"
}

# 显示运行参数
echo ""
echo "运行参数:"
echo "  FASTA文件: $FASTA_FILE"
echo "  PDB目录: $PDB_DIR"
echo "  输出目录: $OUTPUT_DIR"
echo "  生物类型: $ORGANISM_TYPE"
echo "  启用模型: $MODELS"
echo "  使用GPU: $USE_GPU"
echo "  日志级别: $LOG_LEVEL"
echo ""

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 构建命令
CMD="python $PROJECT_DIR/main.py"
CMD="$CMD --fasta $FASTA_FILE"
CMD="$CMD --pdb_dir $PDB_DIR"
CMD="$CMD --output_dir $OUTPUT_DIR"
CMD="$CMD --organism_type $ORGANISM_TYPE"
CMD="$CMD --models $MODELS"
CMD="$CMD --log_level $LOG_LEVEL"

if [ "$USE_GPU" = false ]; then
    CMD="$CMD --cpu_only"
fi

# 运行Pipeline
echo "开始运行Pipeline..."
echo "========================================"
echo ""

eval $CMD

exit_code=$?

echo ""
echo "========================================"
if [ $exit_code -eq 0 ]; then
    echo -e "${GREEN}Pipeline运行成功!${NC}"
    echo ""
    echo "结果文件:"
    ls -la "$OUTPUT_DIR"
else
    echo -e "${RED}Pipeline运行失败，退出码: $exit_code${NC}"
fi

exit $exit_code
