#!/bin/bash

cat $1 |awk '{print $1}' >flist.txt
mkdir ../audio/

for i in `cat flist.txt`;do cp "/home/songhongwei/data_home/DCASE2017-baseline-system/applications/data/TUT-acoustic-scenes-2017-development/$i" ../audio/;done

rm -rf flist.txt

cat $1 |awk '{print $1}' > WavFiles.txt
