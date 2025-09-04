########## Edited on Oct 29th 2024, author: ZX.Yan ##########
/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium
This project is modernized version for event selection,
which can be used to select nuclei(isotopes) with charge 2~8 (From helium to oxygen)

Layout:
include/: Header files declaring functions for events selections
src/: Source files defining functions for events selections
external_libs/: External libraries for RICH beta correction and ROOT dictionaries
main.cpp: Core program for event selection, saving selected events info,
         calculating exposure time, and filling efficiency histograms
build/: CMake build directory and job submission



Compile Guide:

cd /afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/

A. If CMakeLists.txt unchanged:
1. Start bash shell:
   bash
2. Set up environment:
   source ../../useful/amsvar_all.sh
3. Build project:
   cd build
   make -j4

B. If CMakeLists.txt modified:
1. Start bash shell:
   bash
2. Set up environment:
   source ../../useful/amsvar_all.sh
3. Clean and recreate build:
   rm -rf build
   mkdir build
   cd build
4. Configure and build:
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make -j4

cd ../../../submit/

Note: Exit bash shell when checking output rootfile

Running the Program:
After successful compilation, use main_exe:

1. For ISS data:
   ./main_exe <outDir> <outName> <inData> <charge> <UseMass>
   Example:
   ./main_exe yourStorePath/ test.root /eos/ams/group/mit/amsd69n_TMTFNTotHB_B1236P8/1619089051_10.root 3 6 

2. For MC data:
   ./main_exe <outDir> <outName> <inData> "charge|MC" UseMass
   Example:
   ./main_exe yourStorePath/ test.root /eos/ams/group/mit/amsd69n_Li6MCB1306l1_3_6000605N0/1686572662_5.root "3|MC" 6

Parameters:
- outDir: Output directory path
- outName: Output file name
- inData: Input data file path
- charge: Nuclear charge (2-8)
- UseMass: Mass number for binning and converting

/afs/cern.ch/work/s/selu/public/AMS_20230525 find offline here!!

nohup hadd -j 5 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Li6.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/Li6/*.root > & hadd_Li6.log &
nohup hadd -j 5 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Li7.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/Li7/*.root > & hadd_Li7.log &
nohup hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Lit.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Lit/*.root > & hadd_Lit.log &

nohup hadd -j 2 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be7_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/Be7/*.root > & hadd_Be7.log &
nohup hadd -j 4 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be9_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/Be9/*.root > & hadd_Be9.log &
nohup hadd -j 4 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Be10_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/Be10/*.root > & hadd_Be10.log &
nohup hadd -j 5 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B10_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/B10/*.root > & hadd_B10.log &
nohup hadd -j 5 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/B11/*.root > & hadd_B11.log &
nohup hadd -j 4 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/C12_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/C12/*.root > & hadd_C12.log &

nohup hadd -j 6 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7_1.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber7/1{3,4,5}*.root > & shandtxt/hadd_Ber7_1.log &
nohup hadd -j 6 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7_2.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber7/1{6,7}*.root > & shandtxt/hadd_Ber7_2.log &
nohup hadd -j 6 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber9.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber9/*.root > & hadd_Ber9.log &
nohup hadd -j 6 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber10.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber10/*.root > & hadd_Ber10.log &

[zuhao@lxplus993 ~]$ nohup hadd -j 4 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber7.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber7/*.root > & hadd_Ber7.log &
[1] 3975396
[zuhao@lxplus993 ~]$ nohup hadd -j 4 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber9.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber9/*.root > & hadd_Ber9.log &
[2] 3975519
[zuhao@lxplus993 ~]$ nohup hadd -j 4 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Ber10.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber10/*.root > & hadd_Ber10.log &
[3] 3975628
[zuhao@lxplus993 ~]$ ls

nohup hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B10_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/B10/*.root > & hadd_B10.log &
nohup hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B11_bkg.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/B11/*.root > & hadd_B11.log &
nohup hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Bor/*.root > & hadd_Bor.log &


nohup hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/Data_NoAGLCorr/Bor.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Bor/*.root > & /eos/ams/user/z/zuhao/yanzx/Isotope/Data_NoAGLCorr/hadd_Bor.log &
nohup hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/Data_NoAGLCorr/Lit.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Lit/*.root > & /eos/ams/user/z/zuhao/yanzx/Isotope/Data_NoAGLCorr/hadd_Lit.log &
nohup hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/Data_AGLCorr/Lit.root /eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Lit/*.root > & /eos/ams/user/z/zuhao/yanzx/Isotope/Data_AGLCorr/hadd_Lit.log &

hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/Data_NoAGLCorr/C12.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/C12/*.root
hadd -j 10 /eos/ams/user/z/zuhao/yanzx/Isotope/Data_AGLCorr/Li7.root /eos/ams/user/z/zuhao/yanzx/Isotope/MC/Li7/*.root

grep -i "Error" /afs/cern.ch/work/z/zuhao/public/yanzx/logMC/* > check2.log
grep -i "error" /afs/cern.ch/work/z/zuhao/public/yanzx/log9/out* | grep -v "Error in <TChain::SetBranchAddress>: unknown branch" | grep -v "Error in <TTree::SetBranchStatus>: unknown branch"

rm -rf /eos/ams/user/z/zuhao/yanzx/Isotope/MC/{B10,B11,Be7,Be9,Be10,Li6,Li7}/*
