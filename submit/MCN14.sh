#!/bin/bash

# 保存传入的参数
args=("$@")
# 添加库路径
# Source the necessary scripts
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/BashFunction.sh
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/amsvar_all.sh
export LD_LIBRARY_PATH=/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/lib:$LD_LIBRARY_PATH

# 恢复传入的参数
set -- "${args[@]}"

# 输出调试信息，确认接收到的参数
echo "Arguments after sourcing amsvar_all.sh: $@"
echo "Number of arguments: $#"

# Check if a file path is provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <path_to_file> ..."
  exit 1
fi

# Process each file provided as an argument
for line in "$@"; do
    echo "Processing file: $line"
    basename_of_file=$(basename "$line")
    /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/bin/main_exe \
        /eos/user/z/zixuan/Isotope/MC/N14 "$basename_of_file" "$line" "7|MC" 14
done