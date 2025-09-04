#!/bin/bash

# 输入目录
INPUT_DIR="/eos/ams/group/mit/amsd69n_TMTFNTotHB_B1236P8"

# 输出目录
BER7_DIR="/eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber7"
BER9_DIR="/eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber9"
BER10_DIR="/eos/ams/user/z/zuhao/yanzx/Isotope/Data/ISS/Ber10"

# 使用basename来只比较文件名
for dir in "$BER7_DIR" "$BER9_DIR" "$BER10_DIR"; do
    echo "Checking $(basename $dir)..."
    for f in "$INPUT_DIR"/*.root; do
        basename="$(basename "$f")"
        if [ ! -f "$dir/$basename" ]; then
            echo "Missing: $basename"
        fi
    done
done