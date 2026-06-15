#!/bin/bash
# BacterialAntigenFinder Docker构建脚本
#
# 用法:
#   ./build.sh              # 轻量构建（不含模型，运行时挂载）
#   ./build.sh --full       # 完整构建（将模型打包进镜像）
#   ./build.sh --full --no-cache  # 完整构建，不使用缓存
#   ./build.sh --help       # 显示帮助
#
# 说明:
#   轻量构建: 仅包含conda环境和源代码，模型通过运行时卷挂载提供
#            适合开发调试，镜像体积较小
#
#   完整构建: 将模型文件打包进Docker镜像
#            适合部署分发，镜像体积较大
#            需要模型目录与 BacterialAntigenFinder 同级

set -e

# 获取脚本所在目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# 默认参数
BUILD_TYPE="light"
NO_CACHE=""
IMAGE_NAME="bacterial-antigen-finder"
IMAGE_TAG="latest"
DOCKERFILE="${SCRIPT_DIR}/Dockerfile"

# 显示帮助
show_help() {
    echo "BacterialAntigenFinder Docker构建脚本"
    echo ""
    echo "用法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  --full          完整构建（将模型打包进镜像）"
    echo "  --light         轻量构建（默认，不含模型）"
    echo "  --no-cache      不使用Docker缓存"
    echo "  --tag TAG       镜像标签 (默认: latest)"
    echo "  --name NAME     镜像名称 (默认: bacterial-antigen-finder)"
    echo "  -h, --help      显示帮助"
    echo ""
    echo "示例:"
    echo "  # 轻量构建"
    echo "  $0"
    echo ""
    echo "  # 完整构建（需要模型目录在上级目录）"
    echo "  $0 --full"
    echo ""
    echo "  # 完整构建，指定标签"
    echo "  $0 --full --tag v1.0.0"
    echo ""
    echo "  # 不使用缓存重新构建"
    echo "  $0 --full --no-cache"
    echo ""
    echo "目录结构要求（完整构建）:"
    echo "  project_root/"
    echo "  ├── BacterialAntigenFinder/   (本项目)"
    echo "  ├── BepiPred-3.0-main/"
    echo "  ├── DiscoTope-3.0-master/"
    echo "  ├── GraphBepi-main/"
    echo "  └── EpiGraph-main/"
}

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --full)
            BUILD_TYPE="full"
            shift
            ;;
        --light)
            BUILD_TYPE="light"
            shift
            ;;
        --no-cache)
            NO_CACHE="--no-cache"
            shift
            ;;
        --tag)
            IMAGE_TAG="$2"
            shift 2
            ;;
        --name)
            IMAGE_NAME="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "未知参数: $1"
            show_help
            exit 1
            ;;
    esac
done

echo "========================================"
echo "BacterialAntigenFinder Docker构建"
echo "========================================"
echo ""

# 验证项目目录
if [ ! -f "${PROJECT_DIR}/main.py" ]; then
    echo "错误: 未找到项目目录，请确保脚本位于 docker/ 目录下"
    exit 1
fi

# 构建Docker命令
if [ "$BUILD_TYPE" = "full" ]; then
    echo "构建类型: 完整构建（包含模型）"
    echo ""

    # 验证模型目录
    PARENT_DIR="$(dirname "$PROJECT_DIR")"
    MODEL_DIRS=("BepiPred-3.0-main" "DiscoTope-3.0-master" "GraphBepi-main" "EpiGraph-main")

    echo "检查模型目录..."
    MISSING_MODELS=()
    for model in "${MODEL_DIRS[@]}"; do
        if [ -d "${PARENT_DIR}/${model}" ]; then
            echo "  [找到] ${model}"
        else
            echo "  [缺失] ${model}"
            MISSING_MODELS+=("$model")
        fi
    done

    if [ ${#MISSING_MODELS[@]} -gt 0 ]; then
        echo ""
        echo "警告: 以下模型目录缺失: ${MISSING_MODELS[*]}"
        echo "完整构建需要所有模型目录位于: ${PARENT_DIR}/"
        echo ""
        read -p "是否继续构建（缺失模型将不会包含在镜像中）? [y/N] " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "构建已取消"
            exit 0
        fi
    fi

    # 完整构建使用上级目录作为上下文，以便访问模型目录
    BUILD_CONTEXT="${PARENT_DIR}"
    MODEL_SRC_ARG="--build-arg MODEL_SRC=BacterialAntigenFinder/../"

    echo ""
    echo "构建上下文: ${BUILD_CONTEXT}"
    echo "Dockerfile: ${DOCKERFILE}"
    echo "镜像: ${IMAGE_NAME}:${IMAGE_TAG}"
    echo ""

    docker build \
        ${NO_CACHE} \
        -t "${IMAGE_NAME}:${IMAGE_TAG}" \
        -t "${IMAGE_NAME}:full" \
        ${MODEL_SRC_ARG} \
        -f "${DOCKERFILE}" \
        "${BUILD_CONTEXT}"

else
    echo "构建类型: 轻量构建（不含模型）"
    echo ""
    echo "注意: 运行时需要通过 -v 挂载模型目录"
    echo ""
    echo "构建上下文: ${PROJECT_DIR}"
    echo "Dockerfile: ${DOCKERFILE}"
    echo "镜像: ${IMAGE_NAME}:${IMAGE_TAG}"
    echo ""

    docker build \
        ${NO_CACHE} \
        -t "${IMAGE_NAME}:${IMAGE_TAG}" \
        -t "${IMAGE_NAME}:light" \
        -f "${DOCKERFILE}" \
        "${PROJECT_DIR}"
fi

echo ""
echo "========================================"
echo "构建完成!"
echo "镜像: ${IMAGE_NAME}:${IMAGE_TAG}"
echo ""

# 显示使用说明
echo "使用方法:"
echo ""
if [ "$BUILD_TYPE" = "light" ]; then
    echo "  # 运行Pipeline（需挂载模型和数据）"
    echo "  docker run --rm --gpus all \\"
    echo "    -v /path/to/data:/app/data:ro \\"
    echo "    -v /path/to/results:/app/results \\"
    echo "    -v /path/to/BepiPred-3.0-main:/app/models/BepiPred-3.0-main:ro \\"
    echo "    -v /path/to/DiscoTope-3.0-master:/app/models/DiscoTope-3.0-master:ro \\"
    echo "    -v /path/to/GraphBepi-main:/app/models/GraphBepi-main:ro \\"
    echo "    -v /path/to/EpiGraph-main:/app/models/EpiGraph-main:ro \\"
    echo "    ${IMAGE_NAME}:${IMAGE_TAG} \\"
    echo "    --fasta /app/data/antigens.fasta \\"
    echo "    --pdb_dir /app/data/structures/ \\"
    echo "    --organism_type gram- \\"
    echo "    --output_dir /app/results/"
else
    echo "  # 运行Pipeline（模型已内置）"
    echo "  docker run --rm --gpus all \\"
    echo "    -v /path/to/data:/app/data:ro \\"
    echo "    -v /path/to/results:/app/results \\"
    echo "    ${IMAGE_NAME}:${IMAGE_TAG} \\"
    echo "    --fasta /app/data/antigens.fasta \\"
    echo "    --pdb_dir /app/data/structures/ \\"
    echo "    --organism_type gram- \\"
    echo "    --output_dir /app/results/"
fi

echo ""
echo "  # 使用docker-compose运行"
echo "  cd ${SCRIPT_DIR}"
echo "  SPECIES=Klebsiella_pneumoniae docker compose up antigen-finder"
echo ""
echo "  # 运行测试"
echo "  docker compose --profile test up antigen-finder-test"
echo "========================================"
