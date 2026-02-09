#!/bin/bash

# 定义日志和输出文件所在的目录
#LOG_DIR="/afs/cern.ch/work/z/zixuan/logISS/"
LOG_DIR="/afs/cern.ch/work/z/zuhao/public/yanzx/log9/"

# 定义用于存储错误作业列表的输出文件路径 
ERR_JOB_LIST="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist/list_jobs_err2.txt"
# 定义 job list 文件的基础路径
JOB_LIST_BASE_PATH="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist/TxTDivide"

# 初始化错误文件，确保它为空
> "$ERR_JOB_LIST"

# --- 任务 1: 检查 log_ 文件中是否有 "abort" 关键词 ---
echo "=========================================================="
echo "--- 任务 1: 检查 log_* 文件中的 'abort' 关键词 ---"
echo "=========================================================="
FOUND_ABORT=false

for log_file in "$LOG_DIR"/log_*; do
    if [ -f "$log_file" ]; then
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

# --- 任务 2: 检查 output_ 文件的完整性、"successfully" 和 "SysError" ---
echo "=========================================================="
echo "--- 任务 2: 检查 output_* 文件 (需包含 'successfully' 且无 'SysError') ---"
echo "=========================================================="

# 1. 识别所有独特的批次名 (Cluster ID)
BATCH_NAMES=$(find "$LOG_DIR" -maxdepth 1 -type f -name 'output_*.*' | \
               sed 's|.*/output_\([0-9]*\)\..*|\1|' | sort -u)

if [ -z "$BATCH_NAMES" ]; then
    echo "未找到任何 output_*.* 格式的文件。请检查路径或文件名。"
    exit 0
fi

# 2. 对每一个批次名进行检查
EXPECTED_COUNT=2677 # 0 到 2676 共有 2677 个文件

for BATCH_NAME in $BATCH_NAMES; do
    
    echo ""
    echo "----------------------------------------------------------"
    echo ">>> 正在检查批次 (Cluster ID): $BATCH_NAME"
    echo "----------------------------------------------------------"
    
    MISSING_IDS="" 
    SUCCESS_COUNT=0
    FAILURE_COUNT=0
    FAILURE_FILENAMES=""
    RESUBMIT_PIDS="" 
    
    # 检查 0 到 2676 的每个文件
    for i in $(seq 0 2676); do
        FILE_NAME="$LOG_DIR/output_${BATCH_NAME}.${i}"
        
        # A. 检查文件是否存在
        if [ ! -f "$FILE_NAME" ]; then
            MISSING_IDS+="${i} "
            RESUBMIT_PIDS+="${i} "
            continue
        fi
        
        # B. 同时检查 "successfully" (必须有) 和 "SysError" (必须无)
        HAS_SUCCESS=$(grep -q "successfully" "$FILE_NAME"; echo $?)
        HAS_SYSERROR=$(grep -q "SysError" "$FILE_NAME"; echo $?)

        # 逻辑：如果没有 successfully (exit code != 0) OR 有 SysError (exit code == 0)
        if [ $HAS_SUCCESS -eq 0 ] && [ $HAS_SYSERROR -ne 0 ]; then
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        else
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
            # 记录失败原因以便调试
            REASON=""
            [ $HAS_SUCCESS -ne 0 ] && REASON+="[缺少 successfully] "
            [ $HAS_SYSERROR -eq 0 ] && REASON+="[发现 SysError] "
            
            FAILURE_FILENAMES+="$(basename "$FILE_NAME") $REASON\n"
            RESUBMIT_PIDS+="${i} " 
        fi
    done
    
    # 3. 报告结果
    echo "--- 完整性检查报告 ---"
    if [ -z "$MISSING_IDS" ]; then
        echo "文件完整性: 成功。共 $EXPECTED_COUNT 个文件全部存在。"
    else
        echo "文件完整性: 失败。缺少以下 ID: $MISSING_IDS"
    fi

    echo ""
    echo "--- 运行质量检查报告 ---"
    echo "通过 (成功且无系统错误) 的文件数: $SUCCESS_COUNT / $EXPECTED_COUNT"
    echo "失败 (由于错误或未完成) 的文件数: $FAILURE_COUNT / $EXPECTED_COUNT"

    if [ "$FAILURE_COUNT" -gt 0 ]; then
        echo ""
        echo "🚨 警告: 批次 $BATCH_NAME 中有 $FAILURE_COUNT 个作业需要处理。"
        echo -e "$FAILURE_FILENAMES"
    fi
    
    # 4. 写入错误作业列表
    if [ -n "$RESUBMIT_PIDS" ]; then
        echo ""
        echo ">>> 正在记录需要重新提交的作业路径到 $ERR_JOB_LIST"
        for pid in $RESUBMIT_PIDS; do
            ERR_FILE_PATH="${JOB_LIST_BASE_PATH}/job_list_${pid}.txt"
            echo "$ERR_FILE_PATH" >> "$ERR_JOB_LIST"
        done
    fi
done

echo ""
echo "=========================================================="
echo "--- 检查完成 ---"
echo "=========================================================="
