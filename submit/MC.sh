#!/bin/bash

# ================================================================
# 通用 MC 执行脚本 (MC.sh) - V2.0
# 用途：根据 MC 名称和电荷数自动配置参数并执行 SampleProduction
# ================================================================

# ----------------------------------------------------------------
# 1. 读取并存储参数
# ----------------------------------------------------------------
# $1: MCName (e.g., B11, O16)
# $2: Charge (e.g., 5, 8)
# $3: GroupSize (e.g., 80)
# $4: ProcessIndex (e.g., 0, 80, 160, ...)

MCName="$1"
Charge="$2"
GroupSize="$3"
ProcessIndex="$4"

# 检查参数
if [ -z "$MCName" ] || [ -z "$Charge" ] || [ -z "$GroupSize" ] || [ -z "$ProcessIndex" ]; then
  echo "Error: Missing required arguments (MCName, Charge, GroupSize, ProcessIndex)."
  exit 1
fi

echo "--- MC Job Configuration ---"
echo "MC Name: $MCName (Charge: $Charge)"
echo "Group Size: $GroupSize"
echo "Process Index: $ProcessIndex"

# ----------------------------------------------------------------
# 2. 自动配置路径和参数
# ----------------------------------------------------------------

# 2.1 确定输入列表文件路径
INPUT_LIST_DIR="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/submit/joblist"
InputListFile="${INPUT_LIST_DIR}/list_MC_${MCName}.txt"

if [ ! -f "$InputListFile" ]; then
    echo "Error: Input list file not found: $InputListFile"
    exit 1
fi

# 2.2 确定 SampleProduction 目标路径和参数

# 提取质量数 (MassNum)
MassNum=$(echo "$MCName" | sed 's/[^0-9]*//g')

# SampleProduction 参数
PID_PARAM="${Charge}|MC"
EnergyParam="$MassNum" # 假设最后一个参数是质量数

# 目标输出目录 (例如 /eos/user/z/zixuan/Isotope/MC/B10)
OUTPUT_DIR="/eos/user/z/zixuan/Isotope/MC/${MCName}"
# 输出文件名
BaseOutName="MC_${MCName}_Job${ProcessIndex}.root"


echo "Input List File: $InputListFile"
echo "Output Directory: $OUTPUT_DIR"
echo "SampleProduction Params: PID=$PID_PARAM, MassNum=$EnergyParam"

# ----------------------------------------------------------------
# 3. 设置环境
# ----------------------------------------------------------------
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/BashFunction.sh
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/amsvar_all.sh
export LD_LIBRARY_PATH=/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/lib:$LD_LIBRARY_PATH

# ----------------------------------------------------------------
# 4. 文件提取和处理逻辑
# ----------------------------------------------------------------

# 计算起始行和结束行 (Process Index 是起始文件索引)
StartLine=$(( ProcessIndex + 1 ))
EndLine=$(( ProcessIndex + GroupSize ))

echo "Extracting lines from $StartLine to $EndLine"

# 使用 sed 提取指定范围的行
TempFileList=$(mktemp)
sed -n "${StartLine},${EndLine}p" "$InputListFile" > "$TempFileList"

# 检查提取的文件数量
NumFiles=$(wc -l < "$TempFileList")
if [ "$NumFiles" -eq 0 ]; then
    echo "Warning: No files found for this job (Process $ProcessIndex). Exiting."
    rm "$TempFileList"
    exit 0
fi
echo "Successfully extracted $NumFiles files for processing."

# 创建 TChain 输入文件列表
ChainInputFile=$(mktemp --suffix=.txt)
FileCount=0

while IFS= read -r line; do
    if [ ! -z "$line" ]; then
        echo "$line" >> "$ChainInputFile"
        FileCount=$((FileCount + 1))
    fi
done < "$TempFileList"

echo "Created TChain input list: $ChainInputFile containing $FileCount files."

# ----------------------------------------------------------------
# 5. 调用 SampleProduction
# ----------------------------------------------------------------

SAMPLE_PROD_EXE="/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/bin/SampleProduction"

echo "Executing: $SAMPLE_PROD_EXE $OUTPUT_DIR $BaseOutName $ChainInputFile \"$PID_PARAM\" $EnergyParam"

$SAMPLE_PROD_EXE \
    "$OUTPUT_DIR" \
    "$BaseOutName" \
    "$ChainInputFile" \
    "$PID_PARAM" \
    "$EnergyParam"

EXIT_CODE=$?

# ----------------------------------------------------------------
# 6. 清理和退出
# ----------------------------------------------------------------
rm "$TempFileList"
rm "$ChainInputFile"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Job finished successfully. Output: $OUTPUT_DIR/$BaseOutName"
else
    echo "Job FAILED with exit code $EXIT_CODE."
fi

exit $EXIT_CODE