#!/bin/bash

# test the number of input arguments
if [ $# -ne 3 ]
then
	echo "Usage: $0 \"{'WERA_SITE_CODE', ...}\" \"'yyymmddTHH3000'\" \"'yyymmddTHH3000'\""
	exit
fi

date
tic=$(date +%s.%N)
printf "\n"

# Need to set the environment variables relevant for ACORN
source /home/ggalibert/DEFAULT_PATH.env
source /home/ggalibert/STORAGE.env
source /home/ggalibert/ACORN.env

# Run Matlab script to produce FV01 hourly averaged gridded files
matlab -nodisplay -r "cd([getenv('DATA_SERVICES_DIR') '/ACORN']); addpath(fullfile('.', 'Util')); acorn_summary('WERA', true, $1, $2, $3); exit"

printf "\n"
date
toc=$(date +%s.%N)
printf "%6.1Fs\tFV01 radial files processed to hourly averaged gridded files\n"  $(echo "$toc - $tic"|bc )

tic=$(date +%s.%N)
printf "\n"

# Move produced files to OPENDAP
WERA_SOURCE=$DATA/ACORN/WERA/radial_QC/output/datafabric/gridded_1havg_currentmap_QC/
TARGET=$OPENDAP/ACORN/gridded_1h-avg-current-map_QC/
find $WERA_SOURCE -type f -name "*FV01*.nc" -printf %P\\n | sort | rsync -va --remove-source-files --delete-before --files-from=- $WERA_SOURCE $TARGET

printf "\n"
date
toc=$(date +%s.%N)
printf "%6.1Fs\tFV01 hourly averaged gridded files moved from WIP to OPENDAP\n"  $(echo "$toc - $tic"|bc )