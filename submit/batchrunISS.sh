#!/bin/bash

# =========================
# 配置区
# =========================

# 根目录模板（不含 partX）
INPUT_PREFIX="/eos/ams/group/ihep/zixuan/filter/basic_L1Q2p5to8p8_part"
INPUT_SUFFIX=".root"

# part 的范围 [START, END]
START=5
END=9

# 你的 ISSBer.sh 绝对路径
ISSBER="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/ISSBer2.sh"

# time 命令路径
TIME_CMD="/usr/bin/time"

# =========================
# 主流程：只跑一轮
# =========================

echo "==============================="
echo "  Starting single run (no preheat)"
echo "==============================="

for ((i=${START}; i<=${END}; i++)); do
    input_file="${INPUT_PREFIX}${i}${INPUT_SUFFIX}"

    if [[ ! -f "${input_file}" ]]; then
        echo "Warning: input file not found: ${input_file}, skip this one."
        continue
    fi

    log_file="part${i}.log"
    echo "  starting job for ${input_file}, log -> ${log_file}"

    # 后台启动：time 输出和程序输出一起重定向到 log
    (
        "${TIME_CMD}" -v "${ISSBER}" "${input_file}"
    ) &> "${log_file}" &

    echo "    PID: $!"
done

echo
echo "All jobs started in background."
echo "Monitor with:"
echo "  ps -u ${USER} | egrep 'ISSBer.sh|SampleProductio'"
echo
echo "To wait for all jobs in this shell, run:"
echo "  wait"
