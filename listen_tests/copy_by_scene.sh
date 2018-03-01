#!/bin/bash

cat meta.txt |grep $1|awk 'NR<30{print $1}' >flist.txt
mkdir $1

for i in `cat flist.txt`;do cp "/home/songhongwei/data_home/DCASE2017-baseline-system/applications/data/TUT-acoustic-scenes-2017-development/$i" ./$1;done

rm -rf flist.txt

