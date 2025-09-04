#!/bin/bash

src_dir="/eos/ams/group/mit/amsd69n_N15MCB1308l1_7_14000605N0"
output_file="/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/submit/joblist/list_MC_N15.txt"

find "$src_dir" -name "*.root" > "$output_file"

echo "ROOT file paths have been written to $output_file"