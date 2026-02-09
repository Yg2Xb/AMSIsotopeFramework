#!/bin/bash

# =========================
# 参数解析
# =========================

# 检查是否有输入参数
if [ -z "$1" ]; then
    echo "用法: $0 <批次索引(0-3)>"
    echo "  0: 运行 00-04"
    echo "  1: 运行 05-09"
    echo "  2: 运行 10-14"
    echo "  3: 运行 15-19"
    exit 1
fi

BATCH_IDX=$1

# 计算起始和结束编号
START=$((BATCH_IDX * 5))
END=$((START + 4))

# 限制最大范围，防止越界到 part20
if [ $END -gt 19 ]; then
    END=19
fi

# =========================
# 配置区
# =========================

# 你的新 hadd 文件前缀
INPUT_PREFIX="/eos/ams/group/ihep/zixuan/filter/basic_filled_NoHe_part"
INPUT_SUFFIX=".root"

# 你的 ISSBer.sh 绝对路径
ISSBER="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/ISSBer.sh"

# time 命令路径
TIME_CMD="/usr/bin/time"

# =========================
# 主流程
# =========================

echo "==============================="
echo "  Batch Index: ${BATCH_IDX}"
echo "  Range: part$(printf "%02d" ${START}) to part$(printf "%02d" ${END})"
echo "==============================="

for ((i=${START}; i<=${END}; i++)); do
    # 格式化为两位数 (00, 01, ..., 19)
    formatted_i=$(printf "%02d" $i)
    input_file="${INPUT_PREFIX}${formatted_i}${INPUT_SUFFIX}"

    if [[ ! -f "${input_file}" ]]; then
        echo "Warning: [Skip] File not found: ${input_file}"
        continue
    fi

    log_file="part${formatted_i}.log"
    echo "  [Launch] ${input_file} -> ${log_file}"

    # 后台运行
    (
        "${TIME_CMD}" -v "${ISSBER}" "${input_file}"
    ) &> "${log_file}" &

    echo "    PID: $!"
done

echo
echo "Batch ${BATCH_IDX} submitted. Monitor with 'top' or 'ps'."