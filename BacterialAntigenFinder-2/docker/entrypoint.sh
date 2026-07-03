#!/bin/bash
# BacterialAntigenFinder 入口脚本

set -e

# 激活conda
source /opt/conda/etc/profile.d/conda.sh

# 显示物种信息
show_species_info() {
    local species="${SPECIES:-未指定}"
    local organism="${ORGANISM_TYPE:-gram-}"

    echo "物种配置:"
    echo "  目标种类: ${species}"
    echo "  生物类型: ${organism}"

    # 预定义物种信息
    case "${SPECIES}" in
        Pseudomonas_aeruginosa)
            echo "  中文名: 铜绿假单胞菌"
            echo "  耐药性: MDR/XDR"
            ;;
        Klebsiella_pneumoniae)
            echo "  中文名: 肺炎克雷伯菌"
            echo "  耐药性: ESBL/CRE"
            ;;
        Acinetobacter_baumannii)
            echo "  中文名: 鲍曼不动杆菌"
            echo "  耐药性: XDR"
            ;;
        Staphylococcus_aureus)
            echo "  中文名: 金黄色葡萄球菌(MRSA)"
            echo "  耐药性: MRSA"
            ;;
        Escherichia_coli)
            echo "  中文名: 大肠杆菌"
            echo "  耐药性: ESBL"
            ;;
    esac
}

# 验证数据目录
validate_data_dir() {
    if [ ! -d "/app/data" ]; then
        echo "警告: 数据目录 /app/data 不存在，正在创建..."
        mkdir -p /app/data
    fi

    if [ -d "/app/data" ] && [ -z "$(ls -A /app/data 2>/dev/null)" ]; then
        echo "警告: 数据目录 /app/data 为空"
        echo "  请确保已将输入数据挂载到 /app/data 目录"
    fi
}

# 检查模型可用性
check_models() {
    local model_dirs=(
        "/app/models/BepiPred-3.0-main"
        "/app/models/DiscoTope-3.0-master"
        "/app/models/GraphBepi-main"
        "/app/models/EpiGraph-main"
    )

    local available=0
    local missing=0

    echo "模型可用性检查:"
    for dir in "${model_dirs[@]}"; do
        local name=$(basename "$dir")
        if [ -d "$dir" ]; then
            echo "  [可用] $name"
            available=$((available + 1))
        else
            echo "  [缺失] $name"
            missing=$((missing + 1))
        fi
    done

    if [ $missing -gt 0 ]; then
        echo ""
        echo "警告: ${missing} 个模型目录缺失"
        echo "  缺失模型可通过以下方式提供:"
        echo "  1. 运行时通过 -v 挂载模型目录"
        echo "  2. 使用完整构建镜像 (build-arg MODEL_SRC=../)"
        echo "  3. 仅使用可用模型运行 (--models bepipred,discotope)"
    fi

    return $missing
}

# 检查参数
if [ $# -eq 0 ]; then
    echo "BacterialAntigenFinder - 细菌抗原AI智能筛选平台"
    echo ""
    echo "用法:"
    echo "  docker run bacterial-antigen-finder [options]"
    echo ""
    echo "示例:"
    echo "  docker run -v \$(pwd)/data:/app/data -v \$(pwd)/results:/app/results \\"
    echo "      bacterial-antigen-finder \\"
    echo "      --fasta /app/data/antigens.fasta \\"
    echo "      --pdb_dir /app/data/structures/ \\"
    echo "      --organism_type gram- \\"
    echo "      --species Klebsiella_pneumoniae \\"
    echo "      --output_dir /app/results/"
    echo ""
    echo "环境变量:"
    echo "  SPECIES         - 目标细菌种类"
    echo "  ORGANISM_TYPE   - 生物类型 (gram+/gram-/virus)"
    echo ""
    echo "运行 --help 查看完整参数列表"
    exit 0
fi

# 激活主控环境
conda activate master_env

# 显示物种信息
echo "========================================"
show_species_info
echo "========================================"

# 验证数据目录
validate_data_dir

# 检查模型可用性
check_models || true

# 检查环境
echo ""
echo "检查运行环境..."
python -c "import sys; print(f'Python: {sys.version}')"
python -c "import pandas; print(f'Pandas: {pandas.__version__}')"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"

echo ""
echo "启动抗原筛选Pipeline..."
echo "========================================"

# 运行主程序
python /app/main.py "$@"

exit_code=$?

echo "========================================"
if [ $exit_code -eq 0 ]; then
    echo "Pipeline运行完成!"
else
    echo "Pipeline运行失败，退出码: $exit_code"
fi

exit $exit_code
