#!/bin/bash

# 保存传入的参数
args=("$@")

# 添加库路径和环境变量
source /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/useful/BashFunction.sh
source /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/useful/amsvar_all.sh
export LD_LIBRARY_PATH=/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/build/lib:$LD_LIBRARY_PATH

# 恢复传入的参数
set -- "${args[@]}"

# 处理后800行
tail -n +801 /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/submit/joblist/list_MC_O16.txt | while read -r line; do
    echo "Processing: $line"
    /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/build/bin/main_exe \
        /eos/ams/user/z/zuhao/yanzx/Isotope/MC/O16 $(basename "$line") "$line" "8|MC" 16
done
