#!/bin/bash

# This script is an activation script for the python script 20220409_chr_to_whole_genome_merge.py

module load miniconda

DATE=$(date +%H%M_%Y%m%d)
DATA_FOLDER=$1
mkdir "$DATA_FOLDER"/"$DATE"_merge_chr
OUTPUT_FOLDER="$DATA_FOLDER"/"$DATE"_merge_chr

echo "This folder contains run of 20220409_chr_to_whole_genome_merge.py from "$DATE" on the data from "$DATA_FOLDER"" > "$DATA_FOLDER"/"$DATE"_merge_chr/README.md

# submiting:
bsub -q new-short -R "rusage[mem=50000]" -e "$DATA_FOLDER"/"$DATE"_merge_chr/"$DATE"_merge_err.txt -o "$DATA_FOLDER"/"$DATE"_merge_chr/"$DATE"_merge_out.txt conda run -n guy_mmej_env python3 20220409_chr_to_whole_genome_merge.py $DATA_FOLDER $OUTPUT_FOLDER

cp 20220409_chr_to_whole_genome_merge.py "$DATA_FOLDER"/"$DATE"_merge_chr/"$DATE"_20220209_mmej_detection_process.py
