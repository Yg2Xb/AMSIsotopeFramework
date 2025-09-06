#!/bin/bash

# 保存传入的参数
args=("$@")

# 添加库路径和环境变量
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/BashFunction.sh
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/amsvar_all.sh
export LD_LIBRARY_PATH=/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/lib:$LD_LIBRARY_PATH

# 恢复传入的参数
set -- "${args[@]}"

# 处理前800行
head -n 800 /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist/list_MC_O16.txt | while read -r line; do
    echo "Processing: $line"
    /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/bin/main_exe \
        /eos/user/z/zixuan/Isotope/MC/O16 $(basename "$line") "$line" "8|MC" 16
done
