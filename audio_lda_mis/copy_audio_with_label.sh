
cat lda_wrong_table.txt |awk -F ',' 'NR!=1{print $0}' > .full_list.txt
cat lda_wrong_table.txt |awk -F ',' 'NR!=1{print $1".wav"}' > .filelist.txt
cat lda_wrong_table.txt |awk -F ',' 'NR!=1{print $2}' > .truelabel.txt
cat lda_wrong_table.txt |awk -F ',' 'NR!=1{print $3}' > .predlabel.txt


audioroot="/home/songhongwei/data_home/DCASE2017-baseline-system/applications/data/TUT-acoustic-scenes-2017-development/audio/"
for line in `cat .full_list.txt`;do 
	file_name=`echo $line | awk -F ',' '{print $1}'` 	
	truelab=`echo $line |awk -F ',' '{print $2}'`
	predlab=`echo $line |awk -F ',' '{print $3}'`
	if [ "$truelab" == "cafe/restaurant" ];then
		truelab='CafeRestaurant'
	fi

	if [ "$predlab" == "cafe/restaurant" ];then
		predlab='CafeRestaurant'
	fi
	
	cp $audioroot$file_name".wav" ./audio_with_label/$file_name"_"$truelab"_"$predlab".wav";done;
