process=$( qstat | wc -l )
echo $process
while [ $process -eq 5 ] ; do
	sleep 60
	process=$( qstat | wc -l )
done

for var in `seq 1 10` ; do 
	echo -e "\a"
	sleep 0.3
done
