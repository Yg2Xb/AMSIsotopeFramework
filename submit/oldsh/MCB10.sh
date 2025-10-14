#!/bin/bash
#MCB10.sh
# ----------------------------------------------------------------
# 1. 立即读取并存储参数
# ----------------------------------------------------------------
# $1: InputListFile (e.g., /afs/.../list_MC_B10.txt)
# $2: GroupSize (e.g., 10)
# $3: ProcessIndex (e.g., 0, 1, 2, ...)

InputListFile="$1"
GroupSize="$2"
ProcessIndex="$3"

echo "Input List File: $InputListFile"
echo "Group Size: $GroupSize"
echo "Process Index: $ProcessIndex"

# 检查参数
if [ -z "$InputListFile" ] || [ -z "$GroupSize" ] || [ -z "$ProcessIndex" ]; then
  echo "Error: Missing required arguments (InputListFile, GroupSize, ProcessIndex)."
  exit 1
fi

# ----------------------------------------------------------------
# 2. 设置环境（source操作）
# ----------------------------------------------------------------
# 注意：这里不再需要保存和恢复 $args
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/BashFunction.sh
source /afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/useful/amsvar_all.sh
export LD_LIBRARY_PATH=/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/lib:$LD_LIBRARY_PATH

# ----------------------------------------------------------------
# 3. 核心处理逻辑
# ----------------------------------------------------------------

# 计算起始行和结束行
# Condor Process索引从 0 开始
StartLine=$(( (ProcessIndex * GroupSize) + 1 ))
EndLine=$(( StartLine + GroupSize - 1 ))

echo "Extracting lines from $StartLine to $EndLine"

# 使用 sed 提取指定范围的行
TempFileList=$(mktemp)
# 注意：如果 EndLine 超过了文件总行数，sed会停止，不会报错
sed -n "${StartLine},${EndLine}p" "$InputListFile" > "$TempFileList"

# 检查提取的文件数量
NumFiles=$(wc -l < "$TempFileList")
if [ "$NumFiles" -eq 0 ]; then
    echo "Warning: No files found for this job (Process $ProcessIndex). Exiting."
    rm "$TempFileList"
    exit 0
fi
echo "Successfully extracted $NumFiles files for processing."

# 为本次作业创建一个 TChain 输入文件列表
ChainInputFile=$(mktemp --suffix=.txt)
FileCount=0

while IFS= read -r line; do
    if [ ! -z "$line" ]; then
        echo "$line" >> "$ChainInputFile"
        FileCount=$((FileCount + 1))
    fi
done < "$TempFileList"

echo "Created TChain input list: $ChainInputFile containing $FileCount files."

# 确定输出名称
BaseOutName="MC_B10_Job${ProcessIndex}.root"

# 调用 SampleProduction
/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/build/bin/SampleProduction \
    /eos/user/z/zixuan/Isotope/MC/B10 \
    "$BaseOutName" \
    "$ChainInputFile" \
    "5|MC" \
    10

# 清理临时文件
rm "$TempFileList"
rm "$ChainInputFile"

echo "Job finished successfully."