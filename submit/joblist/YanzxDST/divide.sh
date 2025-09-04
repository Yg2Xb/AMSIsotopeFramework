#!/bin/bash

# 定义目录路径
DIR="/eos/ams/user/z/zuhao/yanzx/ams_data/amsd69n_IHEPQcutDST"

# 定义文件计数器和文件组计数器
file_count=0
batch_count=0

# 使用find命令并用sort排序
temp_file=$(mktemp)
find "$DIR" -name "*.root" | sort > "$temp_file"

# 检查是否找到文件
if [ ! -s "$temp_file" ]; then
    echo "No .root files found in the directory."
    rm "$temp_file"
    exit 1
fi

# 逐行读取已排序的文件路径
while IFS= read -r file; do
    # 每10个文件开始一个新的batch
    if (( file_count % 10 == 0 )); then
        # 关闭上一个文件描述符
        if [ $file_count -ne 0 ]; then
            exec 3>&-
        fi
        # 打开新的文件描述符
        exec 3>./job_list_$batch_count.txt
        ((batch_count++))
    fi
    # 写入文件的绝对路径到当前batch文件
    echo "$file" >&3
    ((file_count++))
done < "$temp_file"

# 关闭最后一个文件描述符
exec 3>&-

# 删除临时文件
rm "$temp_file"

echo "Created $batch_count job list files"
echo "Processed $file_count files in total"