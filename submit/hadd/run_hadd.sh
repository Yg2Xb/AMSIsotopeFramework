#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <list_file> <output_name>" >&2
  exit 2
fi

LIST_FILE="$1"
OUT_NAME="$2"

OUTDIR="/eos/user/z/zixuan/Isotope/Add"
LOGDIR="/afs/cern.ch/user/z/zixuan/public/logHadd"
mkdir -p "${OUTDIR}" "${LOGDIR}"

# 控制并行线程，默认 12；可在提交文件中通过 environment 覆盖
THREADS="${HADD_THREADS:-12}"

# 打印环境和调试信息
echo "[run_hadd] Host: $(hostname)"
echo "[run_hadd] Start time: $(date -Iseconds)"
echo "[run_hadd] Threads: ${THREADS}"
echo "[run_hadd] List: ${LIST_FILE}"
echo "[run_hadd] Out : ${OUTDIR}/${OUT_NAME}"

# 打印 ROOT/hadd 版本（若可用）
command -v hadd && hadd -h | head -n 1 || true
command -v root-config && root-config --version || true

# 容错、覆盖、优化、并行
# -f 覆盖，-k 跳过损坏/不可读文件，-O 优化篮子，-j 并行线程
# 使用临时输出避免半成品：写到 .part，再原子 mv
tmp_out="${OUTDIR}/${OUT_NAME}.part.$$"
final_out="${OUTDIR}/${OUTNAME:-${OUT_NAME}}"

# 支持 @list 语法
set -x
hadd -fk -O -j "${THREADS}" "${tmp_out}" @"${LIST_FILE}"
set +x

mv -f "${tmp_out}" "${final_out}"

echo "[run_hadd] Finished at: $(date -Iseconds)"