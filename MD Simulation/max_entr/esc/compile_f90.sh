#!/bin/bash

find . -type f -name "*.f90" | while read file; do
    # 提取文件名（不包括扩展名）
    filename=$(basename "$file" .f90)
    # 提取目录名
    dir=$(dirname "$file")
    # 编译f90文件并将.x文件放在对应的目录
    gfortran "$file" -o "$dir/$filename.x"
done
