#!/bin/bash
rm -r ConstantT
pwd=$(pwd)
find $pwd -type f -name 'f_*' -exec rm -r {} +
#find "`pwd`/Contact_a" -type f -name "*a_[1-9]*" -delete
#find "`pwd`/Contact_a" -type f -name "*a_1[0-9]*" -delete
#find "`pwd`/Contact_a" -type f -name "*a_2[0-9]*" -delete
