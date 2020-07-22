for f in `ls *.txt.gz`; do 
	# echo $f
	#bgzip -f $f &
	tabix -f $f -S 1 -s 3 -b 4 -e 4 &
done
