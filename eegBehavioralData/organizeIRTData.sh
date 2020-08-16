for i in `ls /home/arosen/Documents/bbmcPersonal/eegBehavioralData/careFace/care/*log` ; do
	echo $i
	# Isolate the subject id
	isoString=`echo $i |rev  | cut -f 1 -d / | rev | cut -f 1 -d _`
	# Now cut the final path
	echo $isoString
	# Now create the output files
	outFile1="/home/arosen/Documents/bbmcPersonal/eegBehavioralData/careFace/careMod/pic$isoString.csv"
	outFile2="/home/arosen/Documents/bbmcPersonal/eegBehavioralData/careFace/careMod/resp$isoString.csv"
	## Now grep the files for the rows of interest
	grep "$(printf '\t')Picture" $i | cut -f 1-7  > $outFile1
	grep "$(printf '\t')Response" $i | cut -f 1-7  > $outFile2
done
