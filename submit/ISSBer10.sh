#!/bin/bash

# 保存传入的参数
args=("$@")

# Source the necessary scripts
source /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/useful/BashFunction.sh
source /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/useful/amsvar_all.sh
export LD_LIBRARY_PATH=/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/build/lib:$LD_LIBRARY_PATH

# 恢复传入的参数
set -- "${args[@]}"

# 输出调试信息，确认接收到的参数
echo "Arguments after sourcing scripts: $@"
echo "Number of arguments: $#"

# Check if a file path is provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <path_to_file> ..."
  exit 1
fi

# 逐行读取文件并处理
while IFS= read -r line; do
    echo "Processing file: $line"
     /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/build/bin/main_exe \
        /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber10 \
        "$(basename "$line")" "$line" 4 10
done < "$1"