#!/bin/bash

# ================================================================
# 通用本地批处理脚本 (batchrun.sh) - V2.0
# 用途：计算任务数，控制并发，执行 MC.sh
# ================================================================

# ----------------------------------------------------------------
# 1. 检查输入参数
# ----------------------------------------------------------------
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <MCName> <Charge> <GroupSize>"
    echo "Example: $0 B11 5 80"
    exit 1
fi

MCName="$1"
Charge="$2"
GroupSize="$3"

# ----------------------------------------------------------------
# 2. 配置常量
# ----------------------------------------------------------------
# MC 执行脚本的路径
MC_EXE="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/MC2.sh"
# 输入列表文件路径
INPUT_LIST_FILE="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist/list_MC_${MCName}.txt"
# 节点允许的最大并发数
MAX_PARALLEL_LIMIT=15

echo "--- Batch Run Configuration ---"
echo "Target MC: $MCName (Charge: $Charge)"
echo "Files per Job (GroupSize): $GroupSize"
echo "Max Parallel Processes: $MAX_PARALLEL_LIMIT"
echo "-------------------------------"

# ----------------------------------------------------------------
# 3. 检查文件和计算任务数
# ----------------------------------------------------------------

if [ ! -f "$INPUT_LIST_FILE" ]; then
    echo "Error: Input list file not found: $INPUT_LIST_FILE"
    exit 1
fi

# 计算总行数 (文件总数)
TotalFiles=$(wc -l < "$INPUT_LIST_FILE")
if [ "$TotalFiles" -eq 0 ]; then
    echo "Error: Input list file is empty."
    exit 1
fi

# 计算需要的总作业数
# 使用 bash 算术扩展进行向上取整: ((TotalFiles + GroupSize - 1) / GroupSize)
TotalJobs=$(( (TotalFiles + GroupSize - 1) / GroupSize ))

echo "Total files to process: $TotalFiles"
echo "Total jobs required: $TotalJobs"

# ----------------------------------------------------------------
# 4. 检查并发限制 (仅警告，不停止)
# ----------------------------------------------------------------

if [ "$TotalJobs" -gt "$MAX_PARALLEL_LIMIT" ]; then
    echo "Info: Total jobs ($TotalJobs) is greater than the parallel limit ($MAX_PARALLEL_LIMIT). Using $MAX_PARALLEL_LIMIT concurrency."
fi

# ----------------------------------------------------------------
# 5. 生成 Process Index 列表
# ----------------------------------------------------------------

# 计算最后一个 Process Index
LastIndex=$(( (TotalJobs - 1) * GroupSize ))

# 使用 seq 生成 Process Index 列表
PROCESS_INDICES=$(seq 0 $GroupSize $LastIndex)

# ----------------------------------------------------------------
# 6. 启动后台并发任务
# ----------------------------------------------------------------

LOG_FILE="log_mc_${MCName}_${GroupSize}.2.log"

echo "Starting $TotalJobs jobs with $MAX_PARALLEL_LIMIT concurrency..."
echo "Output log: $LOG_FILE"
echo "-------------------------------"

# 使用 nohup 和 xargs 启动后台任务
# 传递给 MC.sh 的参数顺序: MCName Charge GroupSize ProcessIndex
nohup bash -c "
    echo \"$PROCESS_INDICES\" | xargs -I {} -P $MAX_PARALLEL_LIMIT /usr/bin/time -v $MC_EXE $MCName $Charge $GroupSize {}
" > "$LOG_FILE" 2>&1 &

echo "Batch process started in the background (PID: $!)."
echo "Use 'tail -f $LOG_FILE' to monitor progress."