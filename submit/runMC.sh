#!/bin/bash
# Usage: ./runMC.sh <ParticleName> [mode]
# Example: ./runMC.sh Be7
# mode: "chain" -> chain所有文件一次跑（默认，传入joblist.txt）
#       "loop"  -> 逐文件循环跑

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <ParticleName> [mode]"
  echo "mode: 'chain' or 'loop' (default chain)"
  exit 1
fi

PARTICLE=$1
MODE=${2:-chain}  # 默认 chain
BASE_DIR=/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework
JOBLIST="${BASE_DIR}/submit/joblist/list_MC_${PARTICLE}.txt"
EXEC="${BASE_DIR}/SampleProduction/build/bin/SampleProduction"

# 在 source 前暂时关闭 nounset，防止外部脚本里未定义变量导致退出
was_nounset_set=0
if set -o | grep -q 'nounset[[:space:]]\+on'; then
  was_nounset_set=1
  set +u
fi
# Source 环境
source "${BASE_DIR}/useful/BashFunction.sh"
source "${BASE_DIR}/useful/amsvar_all.sh"
# 恢复 nounset 状态
[ "$was_nounset_set" -eq 1 ] && set -u

# 动态库路径
export LD_LIBRARY_PATH="${BASE_DIR}/SampleProduction/build/lib:${LD_LIBRARY_PATH:-}"

# 输出目录
OUTDIR="/eos/user/z/zixuan/Isotope/MC/${PARTICLE}"
mkdir -p "$OUTDIR"

if [ ! -f "$JOBLIST" ]; then
  echo "Error: Joblist file not found: $JOBLIST"
  exit 1
fi

# 根据粒子名决定 charge (Z) 和质量数 (A)
case "$PARTICLE" in
  B10)  CHARGE=5;  AMASS=10 ;;
  B11)  CHARGE=5;  AMASS=11 ;;
  Be7)  CHARGE=4;  AMASS=7  ;;
  Be9)  CHARGE=4;  AMASS=9  ;;
  Be10) CHARGE=4;  AMASS=10 ;;
  C12)  CHARGE=6;  AMASS=12 ;;
  N14)  CHARGE=7;  AMASS=14 ;;
  O16)  CHARGE=8;  AMASS=16 ;;
  *) echo "Unknown particle: $PARTICLE"; exit 1 ;;
esac

append_root_ext() {
  local name="$1"
  if [[ "$name" == *.root ]]; then
    printf "%s" "$name"
  else
    printf "%s.root" "$name"
  fi
}

echo "Running analysis for $PARTICLE ..."
echo "Output directory: $OUTDIR"
echo "Charge Z=$CHARGE, Mass A=$AMASS"
echo "Mode: $MODE"
echo "Executable: $EXEC"

if [ "$MODE" == "chain" ]; then
    # 传入 joblist.txt，让程序内部展开
    BASENAME="$(append_root_ext "${PARTICLE}_all")"
    echo "Using joblist: $JOBLIST"
    # 简单检查 joblist 中的前几行
    echo "Preview of joblist (first 5 non-comment lines):"
    grep -v '^[[:space:]]*$' "$JOBLIST" | grep -v '^[[:space:]]*#' | head -n 5 || true

    # 执行
    "$EXEC" "$OUTDIR" "$BASENAME" "$JOBLIST" "${CHARGE}|MC" "$AMASS"
else
    # 逐文件循环
    while IFS= read -r line; do
        # 跳过空行和注释
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        base="$(basename "$line")"
        BASENAME="$(append_root_ext "$base")"
        echo "Processing $line -> $BASENAME ..."
        "$EXEC" "$OUTDIR" "$BASENAME" "$line" "${CHARGE}|MC" "$AMASS"
    done < "$JOBLIST"
fi

echo "Analysis for $PARTICLE finished."