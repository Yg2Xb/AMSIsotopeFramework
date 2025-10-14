#!/bin/bash

ZUHAO_PATH="/eos/ams/user/z/zuhao/yanzx/ams_data/amsd69n_IHEPQcutDST"
ZIXUAN_PATH="/eos/user/z/zixuan/Isotope/ISS/Be"
OUTPUT_FILE="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist/Bererr.txt"

# 获取 zuhao 路径下的所有文件名，并创建一个临时文件
find "$ZUHAO_PATH" -type f -printf '%f\n' | sort > /tmp/zuhao_files.txt

# 获取 zixuan 路径下的所有文件名，并创建一个临时文件
find "$ZIXUAN_PATH" -type f -printf '%f\n' | sort > /tmp/zixuan_files.txt

# 使用 comm 命令找出只存在于 zuhao 路径中的文件，即丢失的文件
comm -23 /tmp/zuhao_files.txt /tmp/zixuan_files.txt > /tmp/missing_files.txt

# 遍历丢失的文件名，并找到它们在 zuhao 路径中的绝对路径，然后输出到指定文件
while read -r file; do
    find "$ZUHAO_PATH" -name "$file" >> "$OUTPUT_FILE"
done < /tmp/missing_files.txt

# 清理临时文件
rm /tmp/zuhao_files.txt /tmp/zixuan_files.txt /tmp/missing_files.txt

echo "丢失的文件绝对路径已写入到 $OUTPUT_FILE"