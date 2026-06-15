#!/bin/bash
# BacterialAntigenFinder Docker运行脚本

set -e

echo "========================================"
echo "BacterialAntigenFinder - Docker运行"
echo "========================================"
echo ""

# 获取脚本目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# 默认参数
FASTA_FILE=""
PDB_DIR=""
OUTPUT_DIR="./results"
ORGANISM_TYPE="gram-"
USE_GPU=true
IMAGE_NAME="bacterial-antigen-finder:latest"

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
    echo "  -t, --organism TYPE    生物类型 (默认: gram-)"
    echo "  --cpu                  仅使用CPU"
    echo "  --build                构建Docker镜像"
    echo "  -h, --help             显示帮助"
    echo ""
}

# 构建镜像
build_image() {
    echo "构建Docker镜像..."
    docker build -t "$IMAGE_NAME" -f "$PROJECT_DIR/docker/Dockerfile" "$PROJECT_DIR"
    echo "镜像构建完成!"
}

# 解析参数
BUILD=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--fasta)
            FASTA_FILE="$(realpath "$2")"
            shift 2
            ;;
        -p|--pdb_dir)
            PDB_DIR="$(realpath "$2")"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$(realpath "$2")"
            shift 2
            ;;
        -t|--organism)
            ORGANISM_TYPE="$2"
            shift 2
            ;;
        --cpu)
            USE_GPU=false
            shift
            ;;
        --build)
            BUILD=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "未知参数: $1"
            exit 1
            ;;
    esac
done

# 构建镜像（如果需要）
if [ "$BUILD" = true ]; then
    build_image
fi

# 验证参数
if [ -z "$FASTA_FILE" ] || [ -z "$PDB_DIR" ]; then
    echo "错误: 必须指定FASTA文件和PDB目录"
    show_help
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 构建Docker命令
DOCKER_CMD="docker run --rm"

# GPU支持
if [ "$USE_GPU" = true ]; then
    DOCKER_CMD="$DOCKER_CMD --gpus all"
fi

# 挂载卷
DOCKER_CMD="$DOCKER_CMD -v $FASTA_FILE:/data/input.fasta:ro"
DOCKER_CMD="$DOCKER_CMD -v $PDB_DIR:/data/structures:ro"
DOCKER_CMD="$DOCKER_CMD -v $OUTPUT_DIR:/results"

# 添加模型目录挂载
DOCKER_CMD="$DOCKER_CMD -v $PROJECT_DIR/../BepiPred-3.0-main:/app/models/BepiPred-3.0-main:ro"
DOCKER_CMD="$DOCKER_CMD -v $PROJECT_DIR/../DiscoTope-3.0-master:/app/models/DiscoTope-3.0-master:ro"
DOCKER_CMD="$DOCKER_CMD -v $PROJECT_DIR/../GraphBepi-main:/app/models/GraphBepi-main:ro"
DOCKER_CMD="$DOCKER_CMD -v $PROJECT_DIR/../EpiGraph-main:/app/models/EpiGraph-main:ro"

# 镜像和参数
DOCKER_CMD="$DOCKER_CMD $IMAGE_NAME"
DOCKER_CMD="$DOCKER_CMD --fasta /data/input.fasta"
DOCKER_CMD="$DOCKER_CMD --pdb_dir /data/structures/"
DOCKER_CMD="$DOCKER_CMD --output_dir /results/"
DOCKER_CMD="$DOCKER_CMD --organism_type $ORGANISM_TYPE"

if [ "$USE_GPU" = false ]; then
    DOCKER_CMD="$DOCKER_CMD --cpu_only"
fi

# 显示命令
echo "运行Docker命令:"
echo "$DOCKER_CMD"
echo ""

# 执行
eval $DOCKER_CMD

echo ""
echo "完成! 结果保存在: $OUTPUT_DIR"
