#!/bin/bash

# 20.01.2022

# loading my venv
# module load miniconda
# conda activate guy_mmej_env

# seting the date:
DATE=$(date +%H%M_%Y%m%d)
OUTPUT_FOLDER="/home/labs/alevy/guyta/guy_master_project/results/soybean/Liu_et.al.2020/mmej_detection"

mkdir "$OUTPUT_FOLDER"/"$DATE"
echo "This folder contains the MMEJ detection pipline run from "$DATE"" > "$OUTPUT_FOLDER"/"$DATE"/README.md

# looping the submmition commands per chromosome number
for i in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20}
do
   bsub -q new-medium -R "rusage[mem=50000]" -e "$OUTPUT_FOLDER"/"$DATE"/"$DATE"_mmej_detection_process_err_Chr$i -o "$OUTPUT_FOLDER"/"$DATE"/"$DATE"_mmej_detection_process_output_Chr$i conda run -n guy_mmej_env python3 20220209_mmej_detection_process.py $i $OUTPUT_FOLDER/"$DATE"

   echo "chr$i submitted $DATE"
done

# coping the python script that was used for this run
cp 20220209_mmej_detection_process.py "$OUTPUT_FOLDER"/"$DATE"/"$DATE"_20220209_mmej_detection_process.py