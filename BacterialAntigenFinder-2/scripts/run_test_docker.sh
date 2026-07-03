#!/bin/bash
# BacterialAntigenFinder Docker端到端测试脚本
# 使用Docker容器运行完整Pipeline并验证输出
#
# 用法:
#   ./run_test_docker.sh                     # 运行所有物种测试
#   ./run_test_docker.sh --species pa        # 仅运行铜绿假单胞菌
#   ./run_test_docker.sh --species kp        # 仅运行肺炎克雷伯菌
#   ./run_test_docker.sh --species all       # 运行所有物种(默认)
#   ./run_test_docker.sh --build             # 先构建Docker镜像

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# 获取脚本目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EXAMPLE_DATA_DIR="$PROJECT_DIR/example_data"

# 物种配置
declare -A SPECIES_FASTA
SPECIES_FASTA[pa]="pseudomonas_aeruginosa_antigens.fasta"
SPECIES_FASTA[kp]="klebsiella_pneumoniae_antigens.fasta"

declare -A SPECIES_NAME
SPECIES_NAME[pa]="Pseudomonas_aeruginosa"
SPECIES_NAME[kp]="Klebsiella_pneumoniae"

declare -A SPECIES_DISPLAY
SPECIES_DISPLAY[pa]="铜绿假单胞菌 (Pseudomonas aeruginosa)"
SPECIES_DISPLAY[kp]="肺炎克雷伯菌 (Klebsiella pneumoniae)"

# 默认参数
SPECIES_TO_TEST="all"
OUTPUT_BASE_DIR="$PROJECT_DIR/test_results_docker"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
IMAGE_NAME="bacterial-antigen-finder:latest"
BUILD_IMAGE=false
USE_GPU=true
KEEP_RESULTS=false

# 显示帮助
show_help() {
    echo "BacterialAntigenFinder Docker端到端测试脚本"
    echo ""
    echo "用法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  --species SPEC   指定测试物种 (pa, kp, all, 默认: all)"
    echo "  --output DIR     输出目录 (默认: <project>/test_results_docker)"
    echo "  --build          构建Docker镜像后再运行"
    echo "  --cpu            仅使用CPU (不挂载GPU)"
    echo "  --image NAME     Docker镜像名称 (默认: bacterial-antigen-finder:latest)"
    echo "  --keep           保留测试结果 (默认删除)"
    echo "  -h, --help       显示帮助"
    echo ""
    echo "物种代码:"
    echo "  pa   - 铜绿假单胞菌 (Pseudomonas aeruginosa) [MDR/XDR]"
    echo "  kp   - 肺炎克雷伯菌 (Klebsiella pneumoniae) [ESBL/CRE]"
    echo "  all  - 所有物种 (默认)"
    echo ""
    echo "示例:"
    echo "  $0 --build                      # 构建镜像并运行所有测试"
    echo "  $0 --species pa --cpu           # 仅测试铜绿假单胞菌，不使用GPU"
    echo "  $0 --species kp --keep          # 测试肺炎克雷伯菌并保留结果"
    echo ""
}

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --species)
            SPECIES_TO_TEST="$2"
            shift 2
            ;;
        --output)
            OUTPUT_BASE_DIR="$2"
            shift 2
            ;;
        --build)
            BUILD_IMAGE=true
            shift
            ;;
        --cpu)
            USE_GPU=false
            shift
            ;;
        --image)
            IMAGE_NAME="$2"
            shift 2
            ;;
        --keep)
            KEEP_RESULTS=true
            shift
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

# 确定要测试的物种列表
if [ "$SPECIES_TO_TEST" = "all" ]; then
    TEST_LIST=("pa" "kp")
else
    IFS=',' read -ra TEST_LIST <<< "$SPECIES_TO_TEST"
fi

# 验证物种代码
for sp in "${TEST_LIST[@]}"; do
    if [[ ! -v SPECIES_FASTA[$sp] ]]; then
        echo -e "${RED}错误: 未知物种代码 '$sp'${NC}"
        echo "可用代码: pa, kp"
        exit 1
    fi
done

# 验证测试数据文件存在
for sp in "${TEST_LIST[@]}"; do
    fasta_file="$EXAMPLE_DATA_DIR/${SPECIES_FASTA[$sp]}"
    if [ ! -f "$fasta_file" ]; then
        echo -e "${RED}错误: 测试数据文件不存在: $fasta_file${NC}"
        exit 1
    fi
done

# 检查Docker是否可用
if ! command -v docker &> /dev/null; then
    echo -e "${RED}错误: Docker未安装或不在PATH中${NC}"
    exit 1
fi

# 构建镜像（如果需要）
if $BUILD_IMAGE; then
    echo -e "${BLUE}构建Docker镜像: $IMAGE_NAME${NC}"
    docker build -t "$IMAGE_NAME" -f "$PROJECT_DIR/docker/Dockerfile" "$PROJECT_DIR" || {
        echo -e "${RED}Docker镜像构建失败${NC}"
        exit 1
    }
    echo -e "${GREEN}镜像构建完成!${NC}"
    echo ""
fi

# 检查镜像是否存在
if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo -e "${RED}错误: Docker镜像 '$IMAGE_NAME' 不存在${NC}"
    echo "请使用 --build 参数构建镜像，或指定正确的镜像名称"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_BASE_DIR"

# 测试结果统计
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0
declare -A TEST_RESULTS

echo ""
echo -e "${BLUE}============================================================${NC}"
echo -e "${BLUE}  BacterialAntigenFinder Docker端到端测试${NC}"
echo -e "${BLUE}  时间: $(date '+%Y-%m-%d %H:%M:%S')${NC}"
echo -e "${BLUE}  镜像: $IMAGE_NAME${NC}"
echo -e "${BLUE}  测试物种: ${TEST_LIST[*]}${NC}"
echo -e "${BLUE}============================================================${NC}"
echo ""

# 运行每个物种的测试
for sp in "${TEST_LIST[@]}"; do
    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    species_name="${SPECIES_NAME[$sp]}"
    species_display="${SPECIES_DISPLAY[$sp]}"
    fasta_file="$EXAMPLE_DATA_DIR/${SPECIES_FASTA[$sp]}"
    output_dir="$OUTPUT_BASE_DIR/${sp}_${TIMESTAMP}"

    echo -e "${BLUE}------------------------------------------------------------${NC}"
    echo -e "${BLUE}  测试: $species_display${NC}"
    echo -e "${BLUE}  FASTA: $fasta_file${NC}"
    echo -e "${BLUE}  输出: $output_dir${NC}"
    echo -e "${BLUE}------------------------------------------------------------${NC}"
    echo ""

    # 创建输出目录和PDB输入目录
    mkdir -p "$output_dir"
    pdb_dir="$output_dir/pdb_input"
    mkdir -p "$pdb_dir"

    # 获取绝对路径
    fasta_abs="$(realpath "$fasta_file")"
    output_abs="$(realpath "$output_dir")"
    pdb_abs="$(realpath "$pdb_dir")"
    example_data_abs="$(realpath "$EXAMPLE_DATA_DIR")"

    # 构建Docker命令
    DOCKER_CMD="docker run --rm"

    # GPU支持
    if $USE_GPU; then
        DOCKER_CMD="$DOCKER_CMD --gpus all"
    fi

    # 挂载卷
    DOCKER_CMD="$DOCKER_CMD -v $fasta_abs:/data/input.fasta:ro"
    DOCKER_CMD="$DOCKER_CMD -v $pdb_abs:/data/structures:ro"
    DOCKER_CMD="$DOCKER_CMD -v $output_abs:/results"

    # 挂载模型目录（如果存在）
    models_base="$(dirname "$PROJECT_DIR")"
    for model_dir in BepiPred-3.0-main DiscoTope-3.0-master GraphBepi-main EpiGraph-main Vaxign-ML-docker-master; do
        if [ -d "$models_base/$model_dir" ]; then
            DOCKER_CMD="$DOCKER_CMD -v $models_base/$model_dir:/app/models/$model_dir:ro"
        fi
    done

    # 镜像和Pipeline参数
    DOCKER_CMD="$DOCKER_CMD $IMAGE_NAME"
    DOCKER_CMD="$DOCKER_CMD --fasta /data/input.fasta"
    DOCKER_CMD="$DOCKER_CMD --pdb_dir /data/structures/"
    DOCKER_CMD="$DOCKER_CMD --output_dir /results/"
    DOCKER_CMD="$DOCKER_CMD --organism_type gram-"
    DOCKER_CMD="$DOCKER_CMD --species $species_name"
    DOCKER_CMD="$DOCKER_CMD --output_format both"
    DOCKER_CMD="$DOCKER_CMD --continue_on_error"

    if ! $USE_GPU; then
        DOCKER_CMD="$DOCKER_CMD --cpu_only"
    fi

    # 运行Docker容器
    echo -e "${YELLOW}运行Docker容器...${NC}"
    echo -e "${YELLOW}命令: $DOCKER_CMD${NC}"
    echo ""

    eval $DOCKER_CMD
    exit_code=$?

    if [ $exit_code -ne 0 ]; then
        echo -e "${RED}Docker运行失败 (退出码: $exit_code)${NC}"
        TEST_RESULTS[$sp]="FAILED (exit code: $exit_code)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        continue
    fi

    echo ""
    echo -e "${YELLOW}验证输出文件...${NC}"

    # 验证输出
    validation_passed=true

    # 检查CSV输出
    csv_file="$output_dir/antigen_candidates.csv"
    if [ -f "$csv_file" ]; then
        row_count=$(tail -n +2 "$csv_file" | wc -l)
        echo -e "  ${GREEN}✓${NC} antigen_candidates.csv 存在 ($row_count 行)"
    else
        echo -e "  ${RED}✗${NC} antigen_candidates.csv 不存在"
        validation_passed=false
    fi

    # 检查HTML输出
    html_file="$output_dir/analysis_report.html"
    if [ -f "$html_file" ]; then
        echo -e "  ${GREEN}✓${NC} analysis_report.html 存在"
    else
        echo -e "  ${RED}✗${NC} analysis_report.html 不存在"
        validation_passed=false
    fi

    # 检查元数据JSON
    metadata_file="$output_dir/run_metadata.json"
    if [ -f "$metadata_file" ]; then
        echo -e "  ${GREEN}✓${NC} run_metadata.json 存在"
        # 验证species信息
        if python3 -c "import json; d=json.load(open('$metadata_file')); assert d.get('species') == '$species_name'" 2>/dev/null; then
            echo -e "  ${GREEN}✓${NC} run_metadata.json 包含正确的species信息"
        else
            echo -e "  ${YELLOW}⚠${NC} run_metadata.json species信息不匹配"
        fi
    else
        echo -e "  ${RED}✗${NC} run_metadata.json 不存在"
        validation_passed=false
    fi

    # 记录结果
    if $validation_passed; then
        echo ""
        echo -e "${GREEN}$species_display Docker测试通过!${NC}"
        TEST_RESULTS[$sp]="PASSED"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo ""
        echo -e "${RED}$species_display Docker测试失败!${NC}"
        TEST_RESULTS[$sp]="FAILED (missing output files)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi

    # 清理（除非指定保留）
    if ! $KEEP_RESULTS; then
        echo -e "${YELLOW}清理测试输出...${NC}"
        rm -rf "$output_dir"
    fi

    echo ""
done

# 打印汇总
echo -e "${BLUE}============================================================${NC}"
echo -e "${BLUE}  Docker测试结果汇总${NC}"
echo -e "${BLUE}============================================================${NC}"
echo ""
echo "  镜像: $IMAGE_NAME"
echo "  总测试数: $TOTAL_TESTS"
echo -e "  通过: ${GREEN}$PASSED_TESTS${NC}"
echo -e "  失败: ${RED}$FAILED_TESTS${NC}"
echo ""

for sp in "${TEST_LIST[@]}"; do
    result="${TEST_RESULTS[$sp]}"
    if [[ "$result" == PASSED* ]]; then
        echo -e "  ${GREEN}✓${NC} ${SPECIES_DISPLAY[$sp]}: $result"
    else
        echo -e "  ${RED}✗${NC} ${SPECIES_DISPLAY[$sp]}: $result"
    fi
done

echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}所有Docker测试通过!${NC}"
    exit 0
else
    echo -e "${RED}有 $FAILED_TESTS 个Docker测试失败${NC}"
    exit 1
fi
