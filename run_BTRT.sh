#!/usr/bin/bash

#SBATCH --job-name=BTRT
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --mail-user=david.harrison@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -p swat_plus
#SBATCH -t 48:00:00
#SBATCH --array=0-30

DAYS=($(seq -f "%02g" 1 30))

# cd to directory where job was submitted from
cd $SLURM_SUBMIT_DIR

# get the day information from the array
start_day=${DAYS[$SLURM_ARRAY_TASK_ID]}
end_day=${DAYS[$SLURM_ARRAY_TASK_ID]}

start_month=12
start_year=2014
start_day=31

end_month=$start_month
let end_day=$start_day+1
end_year=$start_year

in_dir=/condo/swatcommon/common/myrorss/$year/$year$start_mon$start_day/segmotion_maxrefl/PolygonTable/scale_0/
out_dir=/condo/swatwork/arkweather/best_track/tracks_RT/201506/segmotion/

# Handle SPC convective days

# 28 day months
if [[ $end_month =~ ^[2]+$ ]]; then
	if [ $end_day = 29 ]; then
		let end_month=$end_month+1
		let end_day=01
	fi			
	#echo "28 days";


# 30 day months
elif [[ $end_month =~ ^[4,6,9,11]+$ ]]; then
	if [ $end_day = 31 ]; then
		let end_month=$end_month+1
		let end_day=01			
	fi
	#echo "30 days";
	
# 31 day months	
elif [[ $end_month =~ ^[1,3,5,7,8,10,12]+$ ]]; then
	if [ $end_day = 32 ]; then
		if [ $end_month = 12 ]; then
			let end_month=01
			let end_day=01
			let end_year=$end_year+1
		else
			let end_month=$end_month+1
			let end_day=01			
		fi
	fi
	#echo "31 days";fi

# Make sure the new month is 0 padded
while [ ${#end_month} -ne 2 ];
	do
		end_month="0"$end_month
	done

echo $SLURM_ARRAY_TASK_ID

#echo $start_month
#echo $start_day
#echo $start_year

#echo $end_month
#echo $end_day
#echo $end_year

python archiveBTRT.py $start_year$start_month$start_day $end_year$end_month$end_day $in_dir $out_dir -it xml -ot seg_json -bt 5 -bd 10 -ht 30






