#!/bin/bash
alias cmake='/cvmfs/sft.cern.ch/lcg/contrib/CMake/latest/Linux-x86_64/bin/cmake'

export wd=`pwd -P`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$wd

export UROOT=6

# Offline directory
export Offline=/cvmfs/ams.cern.ch/Offline
export AMSOPT=/cvmfs/ams.cern.ch/opt

echo "***************************************************"
echo "* Importing amsvar ... "

# icc or gcc
if [[ -z "$UCC" ]]; then
    export UCC=icc64
fi

# Get Redhat version
export USLC=`cat /etc/redhat-release | awk -F'[^0-9]+' '{ print $2 }'`

# Get gcc and icc version with enforcement
if [[ -n "$VGCC" ]]; then
    UGCC=$VGCC
else
    UGCC=0
fi
if [[ -n "$VICC" ]]; then
    UICC=$VICC
else
    UICC=0
fi

# Refine gcc version
if [[ "$UCC" == *gcc* ]]; then
   if [[ "$UGCC" == "0" ]]; then
      if [[ "$USLC" < "7" ]]; then
         UGCC=493
      elif [[ "$USLC" == "7" ]]; then
         UGCC=530
      else
         UGCC=1310
      fi
   fi
else
   if [[ "$UICC" == "0" ]]; then
      if [[ "$USLC" == "5" ]]; then
         UICC=2011
      elif [[ "$USLC" < "8" ]]; then
         UICC=2019
      else
         UICC=2023
      fi
   fi
   if [[ "$UGCC" == "0" ]]; then
      if [[ "$UICC" == "2011" ]]; then
         UGCC=493
      elif [[ "$UICC" == "2019" ]]; then
         UGCC=493
      else
         UGCC=1310
      fi
   fi
fi

# Set compiler versions
export UICC=2023
UGCC=1310

# choose XRDLIB and ROOT
if [[ "$UROOT" == "6" ]]; then
    export XRDLIB=$AMSOPT/xrootd-4.8.2.el9
    export ROOTSYS=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.06/x86_64-almalinux9.4-gcc114-opt
fi

# ROOT and XRDLIB Lib
export LD_LIBRARY_PATH=.:${ROOTSYS}/lib:${XRDLIB}/lib64
export PATH=$ROOTSYS/bin:.:$PATH:/usr/sbin

if [[ "$UROOT" == "6" ]]; then
    pwd=`pwd -P`
    cd $ROOTSYS
    source $ROOTSYS/bin/thisroot.sh
    cd $pwd
fi

# Aux Lib
if [[ "$USLC" > "6" ]]; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/ams.cern.ch/opt/castor.cc7/lib
fi
if [[ "$USLC" == "9" ]]; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/ams.cern.ch/opt/lib64_EL9
fi

# Compiler Lib
export INTEL_LICENSE_FILE=$Offline/intel/licenses

# gcc Lib
if [[ "$UGCC" == "493" ]]; then
    source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh /afs/cern.ch/sw/lcg/external
elif [[ "$UGCC" == "530" ]]; then
    source /cvmfs/sft.cern.ch/lcg/contrib/gcc/5.3.0/x86_64-centos7-gcc53-opt/setup.sh /afs/cern.ch/sw/lcg/external
elif [[ "$UGCC" > "0" ]]; then
    source /cvmfs/sft.cern.ch/lcg/contrib/gcc/11.3.0/x86_64-el9-gcc11-opt/setup.sh /afs/cern.ch/sw/lcg/external
fi

# Intel compiler setup - 修改后的部分
if [[ "$UICC" == "2023" ]]; then
    # 设置 Intel 编译器环境
    export INTELSW=/cvmfs/projects.cern.ch
    export INTELDIR=$INTELSW/intelsw/oneAPI/linux/x86_64/2023
    export INTELVER=compiler/latest/linux
    source $INTELDIR/compiler/latest/env/vars.sh
    
    # 显式设置 Intel 编译器
    export CXX=icpx
    export CC=icx
    export FC=ifx
    
    # 确保 Intel 编译器路径在 PATH 的最前面
    export PATH=$INTELDIR/compiler/latest/linux/bin/intel64:$PATH
elif [[ "$UICC" > "0" ]]; then
    export INTELSW=/cvmfs/projects.cern.ch
    export INTELDIR=$INTELSW/intelsw/oneAPI/linux/x86_64/2024
    export INTELVER=compiler/latest
    source $INTELDIR/$INTELVER/env/vars.sh intel64
fi

# 设置编译标志
export CXXFLAGS="-pthread -std=c++17 -m64 `root-config --cflags` -O3 -D_GLIBCXX_USE_CXX11_ABI=0"
export LDFLAGS="`root-config --ldflags`"

echo "***************************************************"
echo "* USLC=$USLC UCC=$UCC UICC=$UICC UGCC=$UGCC UROOT=$UROOT"
echo "* ROOTSYS=$ROOTSYS"
echo "* Compiler: $CXX"
echo "* CXXFLAGS: $CXXFLAGS"
echo "* `cat /etc/redhat-release`"
echo "***************************************************"

# Verify environment
$CXX --version
echo "ROOT version: $(root-config --version)"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# Add AMSWD environment if needed
if [[ -z "$AMSWD" ]]; then
    export AMSWD=$Offline/vdev
fi

# Add any additional paths specific to your project
if [ -d "/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/external_libs" ]; then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/external_libs
fi