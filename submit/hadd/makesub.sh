#!/usr/bin/env bash
set -euo pipefail

# 基本路径配置
WORKDIR="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/hadd"
LISTDIR="${WORKDIR}/joblist"
LOGDIR="/afs/cern.ch/user/z/zixuan/public/logHadd"
ADD_DIR="/eos/user/z/zixuan/Isotope/Add"

ISS_BE_DIR="/eos/user/z/zixuan/Isotope/ISS/Be"
MC_BASE="/eos/user/z/zixuan/Isotope/MC"
MC_SUBDIRS=("B10" "B11" "Be10" "Be7" "Be9" "C12" "N14" "O16")

# 环境可覆盖
CHUNK_SIZE="${CHUNK_SIZE:-1000}"
HADD_THREADS_DEFAULT="${HADD_THREADS_DEFAULT:-12}"

SUB_FILE="${WORKDIR}/condorHadd.sub"
SUB_FILE_COMPAT="${WORKDIR}/condorHadd_compat.sub"
RUN_SH="${WORKDIR}/run_hadd.sh"

echo "[submit_all] Workdir: ${WORKDIR}"
echo "[submit_all] Listdir: ${LISTDIR}"
echo "[submit_all] Logdir : ${LOGDIR}"
echo "[submit_all] Output : ${ADD_DIR}"
echo "[submit_all] Be dir : ${ISS_BE_DIR}"
echo "[submit_all] MC base: ${MC_BASE}"
echo "[submit_all] CHUNK_SIZE=${CHUNK_SIZE}"
echo "[submit_all] HADD_THREADS_DEFAULT=${HADD_THREADS_DEFAULT}"

mkdir -p "${LISTDIR}" "${LOGDIR}" "${ADD_DIR}"

# 0) 基本检查：run_hadd.sh 可执行
if [[ ! -x "${RUN_SH}" ]]; then
  echo "[INFO] chmod +x ${RUN_SH}"
  chmod +x "${RUN_SH}" || true
fi

# 1) 生成 MC 任务清单
MC_JOBS_FILE="${LISTDIR}/hadd_jobs_mc.txt"
: > "${MC_JOBS_FILE}"
for sub in "${MC_SUBDIRS[@]}"; do
  src="${MC_BASE}/${sub}"
  list="${LISTDIR}/list_${sub}.txt"
  echo "[MC] Listing ${src} -> ${list}"
  find "${src}" -type f -name "*.root" | sort > "${list}" || true
  if [[ ! -s "${list}" ]]; then
    echo "[WARN] ${list} is empty (no ROOT files found)."
  fi
  echo "${list} ${sub}.root" >> "${MC_JOBS_FILE}"
done

# 2) 生成 ISS/Be 分块任务
BE_ALL_LIST="${LISTDIR}/list_Be_all.txt"
echo "[ISS/Be] Listing ${ISS_BE_DIR} -> ${BE_ALL_LIST}"
find "${ISS_BE_DIR}" -type f -name "*.root" | sort > "${BE_ALL_LIST}"
total_files=$(wc -l < "${BE_ALL_LIST}" | awk '{print $1}')
if [[ "${total_files}" -eq 0 ]]; then
  echo "[ERROR] No ROOT files found in ${ISS_BE_DIR}" >&2
  exit 1
fi
echo "[ISS/Be] Total files: ${total_files}"

# 清理旧分块
prefix="${LISTDIR}/list_Be_chunk_"
rm -f ${prefix}*.txt 2>/dev/null || true

# 分块
awk -v n="${CHUNK_SIZE}" -v pfx="${prefix}" '{
  fn=sprintf("%s%03d.txt", pfx, int((NR-1)/n));
  print $0 >> fn
}' "${BE_ALL_LIST}"

# 分块任务清单
BE_JOBS_FILE="${LISTDIR}/hadd_jobs_be_chunks.txt"
: > "${BE_JOBS_FILE}"
for f in ${prefix}*.txt; do
  base=$(basename "$f" .txt)
  chunk_id="${base#list_Be_chunk_}"
  out="Be_chunk_${chunk_id}.root"
  echo "${f} ${out}" >> "${BE_JOBS_FILE}"
done

# 3) 合并总清单
ALL_JOBS_FILE="${LISTDIR}/hadd_jobs_all.txt"
cat "${MC_JOBS_FILE}" "${BE_JOBS_FILE}" > "${ALL_JOBS_FILE}"

echo "[submit_all] Job counts:"
printf "  MC jobs      : %s\n" "$(wc -l < "${MC_JOBS_FILE}" | awk '{print $1}')"
printf "  Be chunk jobs: %s\n" "$(wc -l < "${BE_JOBS_FILE}" | awk '{print $1}')"
printf "  Total jobs   : %s\n" "$(wc -l < "${ALL_JOBS_FILE}" | awk '{print $1}')"

# 4) 校验清单为“两列”
bad_lines=$(awk 'NF!=2 {print NR}' "${ALL_JOBS_FILE}" | wc -l | awk '{print $1}')
if [[ "${bad_lines}" -ne 0 ]]; then
  echo "[ERROR] ${ALL_JOBS_FILE} has lines not having exactly 2 fields." >&2
  awk 'NF!=2 {print NR ": " $0}' "${ALL_JOBS_FILE}" | head -n 10 >&2
  exit 2
fi

# 5) 标准提交文件（不使用文件传输；去掉 when_to_transfer_output）
cat > "${SUB_FILE}" <<'EOF'
Universe  = vanilla
executable  = /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/hadd/run_hadd.sh
arguments   = $(list_file) $(output_name)

output      = /afs/cern.ch/user/z/zixuan/public/logHadd/hadd_out_$(Cluster).$(Process).txt
error       = /afs/cern.ch/user/z/zixuan/public/logHadd/hadd_err_$(Cluster).$(Process).txt
log         = /afs/cern.ch/user/z/zixuan/public/logHadd/hadd_log_$(Cluster).log

+AMSPublic = True
+JobFlavour = "longlunch"

# 不进行文件传输，输入输出都在 AFS/EOS
should_transfer_files   = NO

# 并行线程数由环境变量控制
environment = "HADD_THREADS=12"

queue list_file output_name from /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/hadd/joblist/hadd_jobs_all.txt
EOF

# 6) 兼容版提交文件（旧语法；同样不设置 when_to_transfer_output）
cat > "${SUB_FILE_COMPAT}" <<'EOF'
Universe  = vanilla
executable  = /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/hadd/run_hadd.sh
arguments   = $(item_1) $(item_2)

output      = /afs/cern.ch/user/z/zixuan/public/logHadd/hadd_out_$(Cluster).$(Process).txt
error       = /afs/cern.ch/user/z/zixuan/public/logHadd/hadd_err_$(Cluster).$(Process).txt
log         = /afs/cern.ch/user/z/zixuan/public/logHadd/hadd_log_$(Cluster).log

+AMSPublic = True
+JobFlavour = "longlunch"

should_transfer_files   = NO

environment = "HADD_THREADS=12"

queue 2 from /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/hadd/joblist/hadd_jobs_all.txt
EOF

# 7) 注入 HADD_THREADS 与清单路径
sed -i "s/HADD_THREADS=[0-9][0-9]*/HADD_THREADS=${HADD_THREADS_DEFAULT}/" "${SUB_FILE}"
sed -i "s|${WORKDIR}/joblist/hadd_jobs_all.txt|${ALL_JOBS_FILE}|" "${SUB_FILE}"

sed -i "s/HADD_THREADS=[0-9][0-9]*/HADD_THREADS=${HADD_THREADS_DEFAULT}/" "${SUB_FILE_COMPAT}"
sed -i "s|${WORKDIR}/joblist/hadd_jobs_all.txt|${ALL_JOBS_FILE}|" "${SUB_FILE_COMPAT}"

# 8) 提交
CONDOR_VERSION=$(condor_version 2>/dev/null | head -n1 || echo "unknown")
echo "[submit_all] condor_version: ${CONDOR_VERSION}"
echo "[submit_all] Submitting first wave (MC + Be chunks)..."

if condor_submit "${SUB_FILE}"; then
  echo "[submit_all] Submitted with ${SUB_FILE}"
else
  echo "[submit_all] Fallback to compat submit file. Retrying..."
  condor_submit "${SUB_FILE_COMPAT}"
  echo "[submit_all] Submitted with ${SUB_FILE_COMPAT}"
fi

echo "[submit_all] Done. Monitor with: condor_q $USER"
echo "[submit_all] After jobs finish, to merge Be chunks into one Be.root:"
echo "  ls ${ADD_DIR}/Be_chunk_*.root | sort > ${LISTDIR}/list_Be_chunks_all.txt"
echo "  echo \"${LISTDIR}/list_Be_chunks_all.txt Be.root\" > ${LISTDIR}/hadd_jobs_be_final.txt"
echo "  sed -i \"s|${ALL_JOBS_FILE}|${LISTDIR}/hadd_jobs_be_final.txt|\" ${SUB_FILE}"
echo "  condor_submit ${SUB_FILE}  # 或把 queue 行改成兼容语法: queue 2 from ..."