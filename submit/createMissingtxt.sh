#!/bin/bash

# 定义目录
INPUT_DIR="/eos/ams/group/mit/amsd69n_TMTFNTotHB_B1236P8"
OUTPUT_BASE="/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/submit/joblist"
ISS_BASE="/eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS"

# 确保输出目录存在
mkdir -p "$OUTPUT_BASE"

# 处理Be7、Be9和Be10
for mass in 7 9 10; do
    # 清空或创建新的错误文件
    error_file="$OUTPUT_BASE/Be${mass}err.txt"
    > "$error_file"
    
    # 对比文件并找出缺失的
    while IFS= read -r input_file; do
        basename=$(basename "$input_file")
        if [ ! -f "$ISS_BASE/Ber${mass}/$basename" ]; then
            echo "$INPUT_DIR/$basename" >> "$error_file"
        fi
    done < <(find "$INPUT_DIR" -type f -name "*.root")
    
    # 输出结果统计
    count=$(wc -l < "$error_file")
    echo "Found $count missing files for Be$mass, written to Be${mass}err.txt"
done