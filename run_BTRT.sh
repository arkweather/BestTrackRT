#!/usr/bin/bash

#SBATCH --job-name=BTRT
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --mail-user=david.harrison@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -p swat_plus
#SBATCH -t 48:00:00
#SBATCH --array=0-30

START_DAY=($(seq -f "%02g" 1 30))
END_DAY=($(seq -f "%02g" 1 30))

# cd to directory where job was submitted from
cd $SLURM_SUBMIT_DIR

# get the day information from the array
start_day=${START_DAY[$SLURM_ARRAY_TASK_ID]}
end_day=${END_DAY[$SLURM_ARRAY_TASK_ID]}

end_mon=06
start_mon=06
year=2014
in_dir=/condo/swatcommon/common/myrorss/$year/$year$start_mon$start_day/segmotion_maxrefl/PolygonTable/scale_0/
out_dir=/condo/swatwork/arkweather/best_track/tracks_RT/201506/segmotion/

echo $SLURM_ARRAY_TASK_ID
#echo $in_dir
#echo $(seq -f "%02g" 9 30)

python archiveBTRT.py $year$start_mon$start_day $year$end_mon$end_day $in_dir $out_dir -i xml -t seg_json






