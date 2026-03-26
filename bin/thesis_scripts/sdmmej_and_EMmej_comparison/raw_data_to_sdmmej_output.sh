#!/bin/bash

# ---------------------------------------------------------------------------------
# This sctript is the complete analysis from raw data (*R1.fastq, *R2.fastq)
# to final analysis of the algoritm sdmmej (https://doi.org/10.1093/nar/gkac575)
# ---------------------------------------------------------------------------------

module load BBMap
module load bwa

FULL_DATE=$(date +%H%M_%Y%m%d)
INPUT_R1=$1
INPUT_R2=$2
FILE_RUN_NAME=$3
OURPUT_PATH="${4}/${FULL_DATE}_${3}"
REF=$5

five_prime=$6 # DNA sequence for trim on the 5' end
three_prime=$7 # DNA sequence for trim on the 3' end

mkdir "${OURPUT_PATH}"
ARGS="$@"
echo -e "This folder contains the outputs of: $0\nLocated at: $(pwd) \nTime of submission: $FULL_DATE\nThe command was: $0 $ARGS" > "${OURPUT_PATH}"/"$FULL_DATE"_README.txt

cp /home/labs/alevy/guyta/guy_master_project/scripts/Drosophila/Terrence_Hanscom_NAR_2022/raw_data_to_SAM/raw_data_to_sdmmej_output.sh "${OURPUT_PATH}"
cp /home/labs/alevy/guyta/guy_master_project/scripts/HiFibr/Hi-FiBR-master/Trim_and_Pad.py "${OURPUT_PATH}"
mkdir "${OURPUT_PATH}/bwa_output"
cp /home/labs/alevy/guyta/guy_master_project/scripts/HiFibr/Hi-FiBR-master/Hi-FiBR.py "${OURPUT_PATH}/bwa_output"
cp $REF "${OURPUT_PATH}/bwa_output"

### BBMerge step ###
bsub -q new-short -m public_hosts -R "rusage[mem=10000]" -J "${3}-BBMerge" -e "${OURPUT_PATH}/BBMerge_err.txt" -o "${OURPUT_PATH}/BBMerge_stout.txt" bbmerge.sh in1=${INPUT_R1} in2=${INPUT_R2} out="${OURPUT_PATH}/merged_file/merge.fq" outu1="${OURPUT_PATH}/unmerged_R1.fq" outu2="${OURPUT_PATH}/unmerged_R2.fq" ihist="${OURPUT_PATH}/ihist.txt"
echo "BBMerge step DONE"
sleep 30s

# Triming and padding step
module load miniconda
conda activate guy_mmej_env
# PAD="" # No padding
# echo -e "${OURPUT_PATH}/merged_file\n${five_prime}\n${PAD}\n${three_prime}\n${PAD}" | python3 /home/labs/alevy/guyta/guy_master_project/scripts/HiFibr/Hi-FiBR-master/Trim_and_Pad.py
# -----------------------------------------------------------------------------
# When using Nick's files:
# reference: R0.fa
# cut-suit position from 5' prime: 161 (left)
# cut-suit position from 3' prime: 164 (right)
PAD_5_PRIME="GATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGATCCTAGGAGGGAAAAAATTCGTAC" # TTTGGAGTACGAA
PAD_3_PRIME="TAACGTTAACTCGAGGCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACCCCAGGACC" # GTTAACGTTAACG
echo -e "${OURPUT_PATH}/merged_file\n${five_prime}\n${PAD_5_PRIME}\n${three_prime}\n${PAD_3_PRIME}" | python3 /home/labs/alevy/guyta/guy_master_project/scripts/HiFibr/Hi-FiBR-master/Trim_and_Pad.py
# -----------------------------------------------------------------------------
conda deactivate

### BWA step ###
bwa mem $REF "${OURPUT_PATH}/merged_file/meMatch.fastq" > "${OURPUT_PATH}/bwa_output/bwa_output.sam"
echo "BWA step DONE"

### HiFibr step ###
conda activate guy_mmej_env
# echo -e "${OURPUT_PATH}/bwa_output/R0_WT.fa\n81\n85\n/${OURPUT_PATH}/bwa_output\n30\n30" | python3 "${OURPUT_PATH}/bwa_output/Hi-FiBR.py"
echo -e "${OURPUT_PATH}/bwa_output/R0.fa\n161\n164\n/${OURPUT_PATH}/bwa_output\n30\n30" | python3 "${OURPUT_PATH}/bwa_output/Hi-FiBR.py"

awk -v FS='\t' -v OFS=',' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' "${OURPUT_PATH}/bwa_output/bwa_output_Final.sam" > "${OURPUT_PATH}/bwa_output/bwa_output_Final.csv"

### sdmmej analysis ###

module load miniconda
conda activate sdmmej

HIFIBR_OUTPUT_PATH="${OURPUT_PATH}/bwa_output"

# Adding the header for sdmmej
cat /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Terrence_Hanscom_NAR_2022/from_Nick/Header_for_Hi-FiBR_output_files.csv > "${HIFIBR_OUTPUT_PATH}/bwa_output_Final_with_headers.csv"
echo "\n" >> "${HIFIBR_OUTPUT_PATH}/bwa_output_Final_with_headers.csv"
cat "${HIFIBR_OUTPUT_PATH}/bwa_output_Final.csv" >> "${HIFIBR_OUTPUT_PATH}/bwa_output_Final_with_headers.csv"

cd /home/labs/alevy/guyta/guy_master_project/scripts/sdmmej_algorithm/sdmmej
bsub -q new-short -m public_hosts -R "rusage[mem=10000]" -J sdmmej -e "${HIFIBR_OUTPUT_PATH}/sdmmej_err.txt" -o "${HIFIBR_OUTPUT_PATH}/sdmmej_output.txt" ./run_pipeline.sh "${HIFIBR_OUTPUT_PATH}/bwa_output_Final_with_headers.csv"

# command example:
# When using Nick's files:
# bsub -q new-short -m public_hosts -R "rusage[mem=10000]" -J fastq_to_HiFibr -e "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_WT_data_from_Nick/raw_data_to_sdmmej_output/1023_20220903_err.txt" -o "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_WT_data_from_Nick/raw_data_to_sdmmej_output/1023_20220903_output.txt" ./raw_data_to_sdmmej_output.sh /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Terrence_Hanscom_NAR_2022/from_Nick/R0_WT_R1.fastq /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Terrence_Hanscom_NAR_2022/from_Nick/R0_WT_R2.fastq Iw7  /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_WT_data_from_Nick/raw_data_to_sdmmej_output /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Terrence_Hanscom_NAR_2022/from_Nick/R0.fa TTTGGAGTACGAA GTTAACGTTAACG
# bsub -q new-short -m public_hosts -R "rusage[mem=10000]" -J fastq_to_HiFibr -e "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_POLQ/raw_reads_to_sdmmej_output/1601_20220929_err.txt" -o "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_POLQ/raw_reads_to_sdmmej_output/1601_20220929_output.txt" ./raw_data_to_sdmmej_output.sh /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Terrence_Hanscom_NAR_2022/raw_data/R0_POLQ_reads/SRR19320487_R0_PolQ_1.fastq.gz /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Terrence_Hanscom_NAR_2022/raw_data/R0_POLQ_reads/SRR19320487_R0_PolQ_2.fastq.gz Iw7 /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/R0_POLQ/raw_reads_to_sdmmej_output /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Terrence_Hanscom_NAR_2022/from_Nick/R0.fa TTTGGAGTACGAA GTTAACGTTAACG

