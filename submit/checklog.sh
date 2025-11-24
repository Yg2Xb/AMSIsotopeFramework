#!/bin/bash

# 定义日志和输出文件所在的目录
LOG_DIR="/afs/cern.ch/work/z/zixuan/logISS"

# --- 任务 1: 检查 log_ 文件中是否有 "abort" 关键词 ---
echo "=========================================================="
echo "--- 任务 1: 检查 log_* 文件中的 'abort' 关键词 ---"
echo "=========================================================="
FOUND_ABORT=false

# 查找所有以 log_ 开头的文件
for log_file in "$LOG_DIR"/log_*; do
    # 确保文件存在且是一个常规文件
    if [ -f "$log_file" ]; then
        # 使用 grep 查找 "abort" 并输出匹配行，-q 静默模式
        if grep -q "abort" "$log_file"; then
            echo "在文件 $(basename "$log_file") 中发现 'abort':"
            grep "abort" "$log_file"
            echo "---"
            FOUND_ABORT=true
        fi
    fi
done

if [ "$FOUND_ABORT" = false ]; then
    echo "在所有 log_* 文件中均未发现 'abort' 关键词。"
fi

echo ""
echo ""

# --- 任务 2: 检查 output_ 文件的完整性和 "successfully" 关键词 (针对所有批次) ---
echo "=========================================================="
echo "--- 任务 2: 检查 output_* 文件的完整性和 'successfully' 关键词 (所有批次) ---"
echo "=========================================================="

# 1. 识别所有独特的批次名 (Batch Name)
# 查找所有 output_*.id 文件，并使用 sed 提取出批次号部分 (例如 820886)，然后去重和排序
BATCH_NAMES=$(find "$LOG_DIR" -maxdepth 1 -type f -name 'output_*.*' | \
               sed 's|.*/output_\([0-9]*\)\..*|\1|' | sort -u)

# 如果没有找到任何 output 文件，则退出任务 2
if [ -z "$BATCH_NAMES" ]; then
    echo "未找到任何 output_*.* 格式的文件。请检查路径或文件名。"
    exit 0
fi

# 2. 对每一个批次名进行检查
EXPECTED_COUNT=2676 # 0 到 2950 共有 2951 个文件

for BATCH_NAME in $BATCH_NAMES; do
    
    echo ""
    echo "----------------------------------------------------------"
    echo ">>> 正在检查批次 (Batch Name): $BATCH_NAME"
    echo "----------------------------------------------------------"
    
    MISSING_FILES=""
    SUCCESS_COUNT=0
    FAILURE_COUNT=0
    
    # 检查 0 到 2950 的每个文件
    for i in $(seq 0 2676); do
        FILE_NAME="$LOG_DIR/output_${BATCH_NAME}.${i}"
        
        # 检查文件是否存在
        if [ ! -f "$FILE_NAME" ]; then
            MISSING_FILES+="${i} "
            continue
        fi
        
        # 检查文件中是否包含 "successfully"
        if grep -q "successfully" "$FILE_NAME"; then
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        else
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
            # 可以在这里打印更详细的警告
        fi
    done
    
    # 3. 报告当前批次的检查结果
    
    echo "--- 完整性检查报告 ---"
    
    if [ -z "$MISSING_FILES" ]; then
        echo "文件完整性: 成功。从 0 到 2676 (共 $EXPECTED_COUNT 个文件) 全部存在且连续。"
    else
        echo "文件完整性: 失败。缺少以下 ID 的文件: $MISSING_FILES"
    fi

    echo ""
    echo "--- 'successfully' 关键词检查报告 ---"
    echo "成功 (包含 'successfully') 的文件数: $SUCCESS_COUNT / $EXPECTED_COUNT"
    echo "失败 (缺少 'successfully') 的文件数: $FAILURE_COUNT / $EXPECTED_COUNT"

    if [ "$FAILURE_COUNT" -gt 0 ]; then
        echo "警告: 批次 $BATCH_NAME 中有 $FAILURE_COUNT 个文件没有包含 'successfully'。"
    fi
    
done

echo ""
echo "=========================================================="
echo "--- 所有批次检查完成 ---"
echo "=========================================================="