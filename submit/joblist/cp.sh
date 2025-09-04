#!/bin/bash

# 定义目录路径
dir_bor="/eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber"
dir_lit="/eos/ams/group/mit/amsd69n_TMTFNTotHB_B1236P8"
output_file="/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/submit/joblist/Bererr.txt"

# 清空输出文件
> "$output_file"

# 统计 .root 文件数量
count_bor=$(ls "$dir_bor"/*.root 2>/dev/null | wc -l)
count_lit=$(ls "$dir_lit"/*.root 2>/dev/null | wc -l)

# 输出文件数量到控制台
echo "Ber 文件夹中有 $count_bor 个 .root 文件"
echo "Lit 文件夹中有 $count_lit 个 .root 文件"

# 获取两个目录中的 .root 文件名
files_bor=$(ls "$dir_bor"/*.root 2>/dev/null)
files_lit=$(ls "$dir_lit"/*.root 2>/dev/null)

# 对比文件名并输出到文件
for file in $files_lit; do
    basename_file=$(basename "$file")
    if [[ ! -e "$dir_bor/$basename_file" ]]; then
        echo "$file" >> "$output_file"
    fi
done

# 输出对比结果
echo "对比完成，结果已写入 $output_file"