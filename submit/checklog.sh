#!/bin/bash

# 定义日志和输出文件所在的目录
#LOG_DIR="/afs/cern.ch/work/z/zuhao/public/yanzx/log9/"
LOG_DIR="/afs/cern.ch/work/z/zixuan/logISS/"

# 定义用于存储错误作业列表的输出文件路径 
ERR_JOB_LIST="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist/list_jobs_err2.txt"
# 定义 job list 文件的基础路径 (该路径下的文件命名格式应为 job_list_{ProcessID}.txt)
JOB_LIST_BASE_PATH="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist/TxTDivide"

# 初始化错误文件，确保它为空
> "$ERR_JOB_LIST"

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

# 1. 识别所有独特的批次名 (Cluster ID)
BATCH_NAMES=$(find "$LOG_DIR" -maxdepth 1 -type f -name 'output_*.*' | \
               sed 's|.*/output_\([0-9]*\)\..*|\1|' | sort -u)

# 如果没有找到任何 output 文件，则退出任务 2
if [ -z "$BATCH_NAMES" ]; then
    echo "未找到任何 output_*.* 格式的文件。请检查路径或文件名。"
    exit 0
fi

# 2. 对每一个批次名进行检查
EXPECTED_COUNT=2869 # 0 到 2869 共有 2870 个文件

for BATCH_NAME in $BATCH_NAMES; do
    
    echo ""
    echo "----------------------------------------------------------"
    echo ">>> 正在检查批次 (Cluster ID): $BATCH_NAME"
    echo "----------------------------------------------------------"
    
    MISSING_IDS="" # 存储缺失的 Process ID
    SUCCESS_COUNT=0
    FAILURE_COUNT=0
    FAILURE_FILENAMES=""
    
    # 存储所有需要重新提交的 Process ID (无论是缺失还是失败)
    RESUBMIT_PIDS="" 
    
    # 检查 0 到 2869 的每个文件
    for i in $(seq 0 2869); do
        FILE_NAME="$LOG_DIR/output_${BATCH_NAME}.${i}"
        
        # 检查文件是否存在
        if [ ! -f "$FILE_NAME" ]; then
            MISSING_IDS+="${i} "
            RESUBMIT_PIDS+="${i} " # 缺失的文件 ID 需要重新提交
            continue
        fi
        
        # 检查文件中是否包含 "successfully"
        if grep -q "successfully" "$FILE_NAME"; then
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        else
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
            FAILURE_FILENAMES+="${FILE_NAME}\n"
            RESUBMIT_PIDS+="${i} " # 失败的文件 ID 需要重新提交
        fi
    done
    
    # 3. 报告当前批次的检查结果
    
    echo "--- 完整性检查报告 ---"
    
    if [ -z "$MISSING_IDS" ]; then
        echo "文件完整性: 成功。从 0 到 2869 (共 $EXPECTED_COUNT 个文件) 全部存在且连续。"
    else
        echo "文件完整性: 失败。缺少以下 ID 的文件: $MISSING_IDS"
    fi

    echo ""
    echo "--- 'successfully' 关键词检查报告 ---"
    echo "成功 (包含 'successfully') 的文件数: $SUCCESS_COUNT / $EXPECTED_COUNT"
    echo "失败 (缺少 'successfully') 的文件数: $FAILURE_COUNT / $EXPECTED_COUNT"

    if [ "$FAILURE_COUNT" -gt 0 ]; then
        echo ""
        echo "🚨 警告: 批次 $BATCH_NAME 中有 $FAILURE_COUNT 个文件没有包含 'successfully'。"
        echo "以下是失败的文件列表 (可能需要检查这些文件的内容):"
        echo -e "$FAILURE_FILENAMES"
    fi
    
    # --- 4. 核心功能: 写入错误作业列表 ---
    
    if [ -n "$RESUBMIT_PIDS" ]; then
        
        echo ""
        echo ">>> 正在写入需要重新提交的作业列表..."

        # 遍历所有需要重新提交的 Process ID
        for pid in $RESUBMIT_PIDS; do
            # 构造 job list 文件路径: job_list_{ProcessID}.txt
            ERR_FILE_PATH="${JOB_LIST_BASE_PATH}/job_list_${pid}.txt"
            
            # 将路径写入主错误列表
            echo "$ERR_FILE_PATH" >> "$ERR_JOB_LIST"
            
            echo "记录错误 Job List: $ERR_FILE_PATH"
        done
        
    fi
    
done

echo ""
echo "=========================================================="
echo "--- 所有批次检查完成 ---"
echo "--- 错误 job list 已写入: $ERR_JOB_LIST ---"
echo "=========================================================="
