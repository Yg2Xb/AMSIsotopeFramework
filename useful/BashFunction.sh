
function GetDirName { #Input: TotalFile, Output: Dir, Name

File=$1

File_subdir=` dirname "${File}" `
File_subname=` basename "${File}" `

echo ${File_subdir} ${File_subname}
}

function AMSDisPlay { ###1 InputRoot
currentdir=` pwd `
File=$1
echo "Raw File is: ${File}"
FileList=( $( GetDirName "${File}" ) )
echo "Dir: ${FileList[0]}, Name: ${FileList[1]}" 

if [ ${FileList[0]} == "." ]
then
    File="${currentdir}/${FileList[1]}"
fi
echo "Corrected File is: ${File}"

cd /cvmfs/ams.cern.ch/Offline/vdev/display
source /afs/cern.ch/user/c/czhang/useful/amsvar_icc6_update2017June7AllCVMFS.sh
/afs/cern.ch/ams/Offline/B1100_patches/exe/linuxx8664icc5.34/amsedcPG ${File}

cd "$currentdir"
}

function mkdir_check {
InPut=("$@")
for (( i=0; i < $#; i++ ))
do
if test ! -d ${InPut[$i]}
then
    mkdir -p ${InPut[$i]}
    echo "${InPut[$i]}"
fi
done
}

function CountLengthMax {
#Input: TotalNum NPieces
TotalNum="$1"
NPieces="$2"

NLength=$(( $TotalNum / $NPieces ))

if [ "$(( $TotalNum % $NPieces ))" != 0 ]
then
    (( NLength++ ))
fi

echo $NLength
}

function AutoSetCountLengthMax {
#Input: TotalNum InitalLength MaxNPieces
#Output: NLength NPieces
TotalNum="$1"
InitalLength="$2"
NPiecesMax="$3"

NLength=$InitalLength
NPieces=$( CountLengthMax $TotalNum $NLength )

if [ $NPieces -gt $NPiecesMax ]
then
    NPieces=${NPiecesMax}
    NLength=$( CountLengthMax $TotalNum $NPieces )
fi

echo $NLength $NPieces
}
function CountPieceStartEnd {
#Input: TotalNum NPieces iPiece NLength/""
#OutPut: StartNum EndNum TrueLength
TotalNum=$1
NPieces=$2
iloop=$3

if [ $# -lt 4 ]
then
    NLength=$( CountLengthMax $1 $2 )
else
    NLength=$4
fi

RunStart=$(( ${iloop} * ${NLength} ))
RunEnd=$(( ${iloop} * ${NLength} + ${NLength} - 1 ))

if [ $RunEnd -gt $(( $TotalNum - 1 )) ]
then
    RunEnd=$(( $TotalNum - 1 ))
fi

echo $RunStart $RunEnd $(( $RunEnd - $RunStart + 1 ))
}

function split_setnum { 
#Input: InFile SplitNum OutDir OutName=A

ss_InFile="$1"
ss_Length=` cat ${ss_InFile} | wc -l`
ss_NPieces="$2"
ss_OutDir="$3"
if [ "$4" != "" ]
then
    ss_OutName=$4
else
    ss_OutName="A"
fi

mkdir_check $ss_OutDir

rm -f ${ss_OutDir}/${ss_OutName}*

echo "ss_Length: $ss_Length"
echo "ss_NPieces: $ss_NPieces"

#ss_SplitNum=` expr $ss_Length '/' $ss_NPieces `
#SplitLeft=` expr $ss_Length '%' $ss_NPieces `
ss_SplitNum=$(( $ss_Length / $ss_NPieces ))
SplitLeft=$(( $ss_Length % $ss_NPieces ))
if [ "$SplitLeft" != "0" ]
then
    (( ss_SplitNum++ ))
fi

#split -a 3 -l${ss_SplitNum} $ss_InFile ${ss_OutDir}/${ss_OutName}
split -a 4 -l${ss_SplitNum} -d $ss_InFile ${ss_OutDir}/${ss_OutName} ###ListName is numeric
}

function gvv {
gv "$@" &
}

function FileCheckRun { ### $1: InFileDir || $2: FileReadModel (FILE,*,*1,*List01,...)  || $3: OutDir || $4: OutName || $5: Get Number order 

OutName="ResultCheck.lis"
TmpTXTName="ListofFile.lis"
InFile="$1"
Model="$2"
OutDir="$3"

if [ "$4" != "" ]
then
    OutName="$4"
fi

SepeMark_NumOrder="/ 1"
if [ "$5" != "" ]
then
    SepeMark_NumOrder="$5"
fi

echo "$1 $2 $3"

mkdir_check ${OutDir}

if [ ${Model} == "FILE" ]
then
    echo "InThisSelection"
    cp ${InFile} "${OutDir}/${TmpTXTName}"
    echo "cp ${InFile} ${OutDir}/${TmpTXTName}"
else
    ls ${InFile}${Model}.root > ${OutDir}/${TmpTXTName}
fi

echo "SepeMark_NumOrder is : ${SepeMark_NumOrder}" 
#if [ ${NumOrder} -le 0 ] 
#then
#    NumOrder=1
#fi
####root6 ~/useful/FileCheck.C+'('\"${OutDir}/${TmpTXTName}\"','\"${OutDir}/${OutName}\"','\"${SepeMark_NumOrder}\"')'
root ~/useful/FileCheck.C+'('\"${OutDir}/${TmpTXTName}\"','\"${OutDir}/${OutName}\"','\""${SepeMark_NumOrder}"\"')'
###root ~/useful/FileCheck.C+'('\"${OutDir}/${TmpTXTName}\"','\"${OutDir}/${OutName}\"')'

}
function CombinePiecesRoot { # $1: OutDir, $2: InDir, $3: InName, $4: NCombinePieces

NameEnd="_.root"
if [ "$1" == "" ]
then
    echo "OutDir InDir InName NCombinePieces "
    exit
fi

if [ "$1" != "" ]
then
    OutDir=$1
fi
if [ "$2" != "" ]
then
    InDir=$2
fi
if [ "$3" != "" ]
then
    InName=$3
fi
if [ "$4" != "" ]
then
    NCombinePieces=$4
fi

NAll=` ls ${InDir}/${InName}*${NameEnd} | wc -l ` 

NLength=$( CountLengthMax $NAll $NCombinePieces )

for (( ip=0; ip < $NCombinePieces; ip++ ))
do
    FileList=( "" )
    for(( i=0; i < $NLength; i++ ))
    do
        num=$(( ${NLength} * ${ip} + $i ))
        TmpFile=${InDir}/${InName}${num}${NameEnd}

        if [ -e $TmpFile ]
        then
            FileList=( ${FileList[*]} $TmpFile )
        fi

    done

    hadd -f ${OutDir}/${InName}${ip}${NameEnd} ${FileList[*]}
####    echo "hadd -f ${OutDir}/${InName}${ip}${NameEnd} ${FileList[*]}"

done
}
function GetDataFromListFile { # $1: awk order for data, $2: arraycreatetxt, $3: InFile, $4: OutFile, $5: Options: "recreate"(overwrite outfile), "ForBin"( $1 is the left side of bin, $6 is the right side of bin) 
if [ "$1" == "" ]
then
    echo '$1: awk order for data, $2: arraycreatetxt ( example: double Rbins[] ), $3: InFile, $4: OutFile, $5: Options: "recreate"(overwrite outfile), "ForBin"( $1 is the left side of bin, $6 is the right side of bin)'
    exit
fi

NumOrder="$1"
ArrayTxt="$2"
InFile="$3"
OutFile="$4"
Option="$5"

if [[ "$Option" == *"recreate"* ]]
then
    echo "${ArrayTxt} = { " > ${OutFile}
else
    echo "${ArrayTxt} = { " >> ${OutFile}
fi

Array=( $(cat ${InFile} | awk '{printf " %s ", '\$${NumOrder}' }') )
Length=${#Array[*]}

if [[ "$Option" == *"ForBin"* ]]
then
    NumOrderRight="$6"    
    BinRight=( $(cat ${InFile} | awk '{printf " %s ", '\$${NumOrderRight}' }') )
    LengthRight=${#BinRight[*]}
    EndBin=$(( ${LengthRight} - 1 ))

    Array[$Length]=${BinRight[$EndBin]}
    Length=${#Array[*]}
fi

for (( i=0; i < ${Length}; i++ ))
do
    if [[ "${ArrayTxt}" != "TString"* ]]
    then
        if [ "$i" == "0" ]
        then
            echo -ne "${Array[$i]} " >> ${OutFile}
        else
            echo -ne ",\t ${Array[$i]}" >> ${OutFile}
        fi
    else
        if [ "$i" == "0" ]
        then
            echo -ne "\" ${Array[$i]} \"" >> ${OutFile}
        else
            echo -ne ",\t \" ${Array[$i]} \"" >> ${OutFile}
        fi
    fi
done
echo "}; " >> ${OutFile}
}

function TransFiles { ### $1: File List, $2: Source Dir, $3: Target Dir

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ]
then
    echo "1: File List, 2: Source Dir, 3: Target Dir, 4: Option 'yes-go'---no reconfirm, "
    exit
fi

FileList=$1
SourceDir=$2
TargetDir=$3
Option="$4"

if [ "${Option}" != "yes-go" ]
then
    echo " ==== Confirm ===="
    echo "File List: ${FileList}"
    echo "!!! Please Confirm carefully !!!"
    echo "SourceDir: ${SourceDir}"
    echo "TargetDir(File will be remove and copy from SourceDir): ${TargetDir}"
    echo "y or n: "

    read JudgeWord
    if [ "${JudgeWord}" != "y" ]
    then
        exit
    fi
fi

cat ${FileList} | xargs -I {} basename {} | xargs -I {} rm -f ${TargetDir}/{}
cat ${FileList} | xargs -I {} basename {} | xargs -I {} eos cp ${SourceDir}/{} ${TargetDir}/

}
function lstofile { ### $1: InDir, $2: OutFile, $3: InName/InKind(Inital: "*.root")

if [ "$1" == "" ] || [ "$2" == "" ]
then
    echo "1: InDir, 2: OutFile, 3: InName/InKind(Inital: *.root)"
    exit
fi

ltf_InDir=$1
ltf_OutFile="$2"

ltf_OutDir=$( dirname ${OutFile} )
mkdir_check $ltf_OutDir


ltf_InKind="*.root"
if [ "$3" != "" ]
then
ltf_InKind="$3"
fi

find ${ltf_InDir}"/" -type f -name "${ltf_InKind}" -print0 | xargs -0 ls > ${ltf_OutFile}

}
function haddQuick { ### $1: OutFile, $2: InFiles, $3 Pieces: LimitMax is 10(need be improved)
    
OutFile="$1"
InFiles="$2"
NPiece=$3

currentdir=` pwd `

#rm -f ${OutFile}_tmp*.root
#rm -f ${OutFile}.root

echo "${currentdir}"

for((i=0;i<${NPiece};i++))
do

    echo xterm -T "HAddPiece${i}" -e "hadd -f ${OutFile}_tmp_${i}_.root ${InFiles}*${i}_.root ; sleep 10"
    sleep 5
done | xargs -l1 -P${NPiece} -i sh -c "{}"

sleep 5

hadd -f "${OutFile}.root" ${OutFile}_tmp*.root

rm -f ${OutFile}_tmp*.root

}
function CreateLineJobList { ### $1: InDir, $2: OutDir, $3: OutName, $4: NPieces, $5 Option: "ReadExist"

InDir="$1"
OutDir="$2"
OutName="$3"
NPieces="$4"
Option="$5"

OutFile_AllList="${OutDir}/AllFileList_${OutName}.lis"
OutFile_subdirlist="${OutDir}/submit_${OutName}/"
OutFile_test="${OutDir}/${OutName}_test.lis"

OutFile="${OutDir}/${OutName}.lis"

##Debug##
echo "OutFile_AllList: ${OutFile_AllList}"
echo "OutFile: ${OutFile}"
##Debug##
TestFileNum=$(( ${NPieces} - 1 ))
Flag_ReRead=true

if [[ "${Option}" == *"ReadExist"* ]]
then
    Flag_ReRead=false
fi

if ${Flag_ReRead}
then
    lstofile "${InDir}" "${OutFile_AllList}"
fi

split_setnum "${OutFile_AllList}" "${NPieces}" "${OutFile_subdirlist}" "A"

##Debug##
echo "OutFile_AllList: ${OutFile_AllList}"
echo "OutFile: ${OutFile}"
##Debug##

ls ${OutFile_subdirlist}/A* > ${OutFile}

##Debug##
echo "OutFile_AllList: ${OutFile_AllList}"
echo "OutFile: ${OutFile}"
##Debug##

cat ${OutFile_AllList} | awk ' NR == '${TestFileNum}' {print $0 } ' > ${OutFile_test}
}

function subjob_condor { ### $1: jobaddname(normally is No. ), $2: jobandlogsavedir, $3: dir/exe, $4: exe inputs(need all in ""), $5: environment_need_source, $6: dirTotrans, $7: filetotrans, $8: dirraw
####Attention: if want to tmp save root at service, outdir in $2 should be "\${TMPDIR}", $5 is the finnal dir.

if [ "${1}" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ]
then
    echo " === UseGuid === " 
    echo "Must Fill: 1: jobaddname(normally is No. ), 2: jobandlogsavedir, 3: dir/exe, 4: exe inputs(need all in "") "
    echo "Not Must: 5: environment_need_source"
    echo "Not Must: 6: dirTotrans, 7: filetotrans, 8: dirraw"
    echo " =============== "
    exit
fi

sc_AddName=${1}
sc_JobSubDir=${2}
sc_RunEXE=${3}
sc_RunInputs=${4}

if [ "$5" != "" ]
then
    sc_evn_source="${5} ;"
else
    sc_evn_source=""
fi

sc_FlagTrans=false
if [ "$6" != "" ] && [ "$7" != "" ] && [ "$8" != "" ]
then
    sc_DirToTrans=${6}
    sc_FileToTrans=${7}
    sc_DirRaw=${8}

    sc_FlagTrans=true
fi

chmod +x ${sc_RunEXE}
mkdir_check ${sc_JobSubDir}/joblog
mkdir_check ${sc_DirToTrans}

echo "#!/bin/bash" > ${sc_JobSubDir}/RunJob_${sc_AddName}.sh

sc_RunPart="date; ${sc_evn_source} ${sc_RunEXE} ${sc_RunInputs} ; date"

if ${sc_FlagTrans}
then
echo "${sc_RunPart} ; eos cp ${sc_DirRaw}/${sc_FileToTrans}* ${sc_DirToTrans}/ || eos cp ${sc_DirRaw}/${sc_FileToTrans}* ${sc_DirToTrans}/ " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sh
else
    echo "${sc_RunPart}" >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sh
fi

chmod +x ${sc_JobSubDir}/RunJob_${sc_AddName}.sh

echo "executable = ${sc_JobSubDir}/RunJob_${sc_AddName}.sh " > ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "arguments = \"\$(Cluster) \$(Process) \" " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "output = ${sc_JobSubDir}/joblog/Run_${sc_AddName}_\$(Process).out " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "error = ${sc_JobSubDir}/joblog/Run_${sc_AddName}_\$(Process).err " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "log = ${sc_JobSubDir}/joblog/Run_${sc_AddName}_\$(Process).log " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "transfer_output_files = \"\" " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "transfer_input_files = \"\" " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "+JobFlavour = \"longlunch\" " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "getenv = true " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
echo "queue " >> ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
chmod +x ${sc_JobSubDir}/RunJob_${sc_AddName}.sub

condor_submit ${sc_JobSubDir}/RunJob_${sc_AddName}.sub
}

function subjob_condor_InLine { ### $1: scI_ListFile, others same as subjob_condir

if [ "$1" == "" ]
then
    echo "<====== UseGuide ======>"
    echo "1: scI_ListFile instand of jobaddname"
    echo "others same as subcondir, which is: "
    subjob_condir ""
    echo "<======================>"
fi

scI_Inputs=("$@")
scI_Inputs_used=("$@")

scI_ListFile=${scI_Inputs[0]}
scI_JobSubDir=${scI_Inputs[1]}

rm -rf ${scI_JobSubDir}

scI_icount=0

while read InFile
do
    scI_Inputs_used[0]=${scI_icount}
    echo"Run: ${scI_Inputs_used[*]}"

    subjob_condir ${scI_Inputs_used[*]}

    (( scI_icount++ ))
done < ${scI_ListFile}
}

