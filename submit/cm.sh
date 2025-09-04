#!/bin/bash

base="/afs/cern.ch/work/z/zuhao/public/yanzx/log9/output_7271180."
start=0
end=2460

echo "开始检查缺失文件..."
echo "从 ${base}${start} 到 ${base}${end}"

missing=0
for i in $(seq $start $end); do
    if [ ! -f "${base}${i}" ]; then
        echo "缺失: ${base}${i}"
        ((missing++))
    fi
done

echo "检查完成"
echo "总共缺失文件数: $missing"