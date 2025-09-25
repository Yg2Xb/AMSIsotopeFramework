#!/usr/bin/env bash
set -euo pipefail

# 配置
ISS_BE_DIR="/eos/user/z/zixuan/Isotope/ISS/Be"
MC_BASE="/eos/user/z/zixuan/Isotope/MC"
MC_SUBDIRS=("B10" "B11" "Be10" "Be7" "Be9" "C12" "N14" "O16")
ADD_DIR="/eos/user/z/zixuan/Isotope/Add"

# 工作区与列表输出
WORKDIR="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/hadd"
LISTDIR="${WORKDIR}/joblist"
LOGDIR="/afs/cern.ch/user/z/zixuan/public/logHadd"

mkdir -p "${LISTDIR}" "${LOGDIR}"
mkdir -p "${ADD_DIR}"

echo "[make_hadd_lists] Generating file lists in ${LISTDIR}"

# 1) MC 各子目录列表与任务清单
MC_JOBS_FILE="${LISTDIR}/hadd_jobs_mc.txt"
: > "${MC_JOBS_FILE}"

for sub in "${MC_SUBDIRS[@]}"; do
  src="${MC_BASE}/${sub}"
  list="${LISTDIR}/list_${sub}.txt"
  echo "[MC] ${sub} -> ${list}"
  # 仅收集普通文件，排序可稳定化 IO
  find "${src}" -type f -name "*.root" | sort > "${list}"
  # 输出目标文件名
  echo "${list} ${sub}.root" >> "${MC_JOBS_FILE}"
done

# 2) ISS/Be 分块
BE_ALL_LIST="${LISTDIR}/list_Be_all.txt"
echo "[ISS/Be] Listing all ROOT files -> ${BE_ALL_LIST}"
find "${ISS_BE_DIR}" -type f -name "*.root" | sort > "${BE_ALL_LIST}"

# 分块大小，可按需要调整（1000 左右较稳妥）
CHUNK_SIZE=1000

# 将 BE_ALL_LIST 按行分块
# 生成 list_Be_chunk_000.txt, 001.txt, ...
prefix="${LISTDIR}/list_Be_chunk_"
# 清理历史块
rm -f ${prefix}*.txt 2>/dev/null || true

total_files=$(wc -l < "${BE_ALL_LIST}" | awk '{print $1}')
if [[ "${total_files}" -eq 0 ]]; then
  echo "ERROR: No ROOT files found in ${ISS_BE_DIR}" >&2
  exit 1
fi

# 使用 awk 分块，编号从 0 开始
awk -v n="${CHUNK_SIZE}" -v pfx="${prefix}" '{
  fn=sprintf("%s%03d.txt", pfx, int((NR-1)/n));
  print $0 >> fn
}' "${BE_ALL_LIST}"

# 生成 Be 一级合并任务清单
BE_JOBS_FILE="${LISTDIR}/hadd_jobs_be_chunks.txt"
: > "${BE_JOBS_FILE}"
for f in ${prefix}*.txt; do
  base=$(basename "$f" .txt)           # list_Be_chunk_000
  chunk_id="${base#list_Be_chunk_}"    # 000
  out="Be_chunk_${chunk_id}.root"
  echo "${f} ${out}" >> "${BE_JOBS_FILE}"
done

# 生成 Be 二级合并列表（由 submit_be_final.sh 使用）
BE_CHUNKS_LIST="${LISTDIR}/list_Be_chunks_all.txt"
ls "${ADD_DIR}"/Be_chunk_*.root 2>/dev/null | sort > "${BE_CHUNKS_LIST}" || true

echo "[make_hadd_lists] Done."
echo "MC jobs list: ${MC_JOBS_FILE}"
echo "Be chunk jobs list: ${BE_JOBS_FILE}"
echo "Be chunks existing list (may be empty before running): ${BE_CHUNKS_LIST}"