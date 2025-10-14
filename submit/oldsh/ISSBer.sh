#!/bin/bash

# 保存传入的参数，避免 source 过程中变量干扰
args=("$@")

# 加载环境（不使用 set -u，避免 amsvar_all.sh 中未定义变量导致中断）
# 如需更稳妥，可在此处临时关闭 nounset： set +u; source ...; set -u
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/BashFunction.sh
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/amsvar_all.sh

# 本项目库路径
export LD_LIBRARY_PATH=/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/lib:$LD_LIBRARY_PATH

# 恢复传入的参数
set -- "${args[@]}"

# 调试输出
echo "Arguments after sourcing scripts: $@"
echo "Number of arguments: $#"

# 参数检查：期望传入一个“子列表 txt”
if [ $# -lt 1 ]; then
  echo "Usage: $0 <list_txt>"
  exit 1
fi

LISTFILE="$1"

# 输出目录与文件名（基于子列表名生成，确保一作业一个输出）
OUTDIR="/eos/user/z/zixuan/Isotope/ISS/Be"
BASENAME="$(basename "${LISTFILE%.*}")"   # 例如 job_list_4
OUTNAME="${BASENAME}.root"                # 例如 job_list_4.root

echo "[INFO] Using list: $LISTFILE"
echo "[INFO] Output file: $OUTDIR/$OUTNAME"

# 一次性运行：将列表文件作为 inData 传入，程序内部 TChain 读取多个 ROOT 并写出一个输出
/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/bin/SampleProduction \
  "$OUTDIR" "$OUTNAME" "$LISTFILE" 4 7

echo "[INFO] Finished: $OUTDIR/$OUTNAME"