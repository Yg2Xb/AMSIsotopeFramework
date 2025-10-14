#!/bin/bash

# 定义日志路径
LOG_DIR1="/afs/cern.ch/work/z/zixuan/logMC"
LOG_DIR2="/afs/cern.ch/work/z/zixuan/logISS"

# 输出文件名
ERROR_LOG="check2.log"
MISSING_LOG="missing_files.log"

# 清空旧日志
> "$ERROR_LOG"
> "$MISSING_LOG"

# 功能1: 收集错误日志（修复路径问题 + 扩展忽略规则）
{
  # 从 logMC 下所有 .log 中收集含 Error 的行，并忽略若干已知无害模式
  grep -i "Error" "$LOG_DIR1"/*.log 2>/dev/null \
    | grep -v 'Error in <TChain::SetBranchAddress>: unknown branch' \
    | grep -v 'Error in <TTree::SetBranchStatus>: unknown branch' \
    | grep -v 'Error in <GSLError>' \
    | grep -v 'Error found in integrating function'

  # 从 logISS 下 log_* 中收集含 error 的行，并同样忽略无害模式
  grep -i "error" "$LOG_DIR2"/log_* 2>/dev/null \
    | grep -v 'Error in <TChain::SetBranchAddress>: unknown branch' \
    | grep -v 'Error in <TTree::SetBranchStatus>: unknown branch' \
    | grep -v 'Error in <GSLError>' \
    | grep -v 'Error found in integrating function'
} > "$ERROR_LOG"

# 功能2: 检查不连续的 output_*.i 文件（按任务号分组检查）
check_sequence() {
  local dir=$1
  local prefix=$2

  find "$dir" -name "${prefix}_*.*" -type f | while read -r file; do
    taskid=$(basename "$file" | cut -d_ -f2 | cut -d. -f1)
    i=$(basename "$file" | cut -d. -f2)
    echo "$taskid $i"
  done | sort -k1,1n -k2,2n | uniq | while read -r taskid i; do
    echo "$taskid" >> "${dir}/.taskids.tmp"
  done

  sort -n "${dir}/.taskids.tmp" | uniq | while read -r taskid; do
    is=$(find "$dir" -name "${prefix}_${taskid}.*" -type f \
         | xargs -n1 basename \
         | cut -d. -f2 \
         | sort -n)

    prev=-1
    for i in $is; do
      if [ "$prev" -ne -1 ] && [ "$i" -ne $((prev + 1)) ]; then
        echo "[WARNING] 在目录 $dir 任务 $taskid 中发现不连续i值: 缺失 ${prefix}_${taskid}.$((prev + 1))" >> "$MISSING_LOG"
      fi
      prev=$i
    done

    first_i=$(echo "$is" | head -n1)
    if [ -n "$first_i" ] && [ "$first_i" -ne 0 ]; then
      echo "[WARNING] 在目录 $dir 任务 $taskid 中起始i值不是0: 第一个i是 $first_i" >> "$MISSING_LOG"
    fi
  done

  rm -f "${dir}/.taskids.tmp"
}

# 检查两个目录的output文件
check_sequence "$LOG_DIR1" "output"
check_sequence "$LOG_DIR2" "output"

# 输出结果
echo "===== 错误日志已保存到 $ERROR_LOG ====="
if [ -s "$MISSING_LOG" ]; then
  echo "===== 发现不连续文件 ====="
  cat "$MISSING_LOG"
else
  echo "===== 所有i值连续 ====="
fi