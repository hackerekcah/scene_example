
cat lda_wrong_table.txt |awk -F ',' 'NR!=1{print $1".wav"}' > .filelist.txt
audioroot="/home/songhongwei/data_home/DCASE2017-baseline-system/applications/data/TUT-acoustic-scenes-2017-development/audio/"
for file in `cat .filelist.txt`;do cp $audioroot$file ./audio/;done;
