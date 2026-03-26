#!/bin/bash

# This is a pipline that operates the different stages of processing and analyzing
# datasets in the project

module load bedtools
module load miniconda
conda activate guy_mmej_env

DATE=$(date +%Y%m%d)
FULL_DATE=$(date +%H%M_%Y%m%d)
PATH_TO_OUTPUT_DATA=$1
PATH_TO_DATA=$2
DATA_FILE_NAME=$3
REF_GENOME_FILE=$4
FULL_DATE=$5
CONTEXT_WINDOW_SIZE=$6 # 1 based window size, change this line per run
# SIMS=$7

# crearing a folder for the run's output
mkdir "$PATH_TO_OUTPUT_DATA"/"$FULL_DATE"_main_pipline_output_"$DATA_FILE_NAME"
PATH_TO_OUTPUT_DATA="$PATH_TO_OUTPUT_DATA"/"$FULL_DATE"_main_pipline_output_"$DATA_FILE_NAME"
PATH_TO_REALINGED_OUTPUT_DATA="$PATH_TO_OUTPUT_DATA"

mkdir $PATH_TO_OUTPUT_DATA/scripts
mkdir $PATH_TO_OUTPUT_DATA/docs
# copying all relevant scripts to the output/scripts folder
cp 20220503_main_pipline.sh "$PATH_TO_OUTPUT_DATA"/scripts
cp 20220503_main_pipline_submition.sh "$PATH_TO_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/src/MicroHomology_module_v3.py "$PATH_TO_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/src/MMEJ_2nd_order_MM_v2.py "$PATH_TO_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/src/ExpectationMaximization.py "$PATH_TO_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/src/data_exploration_util.py "$PATH_TO_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/src/pysam_getfasta.py "$PATH_TO_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/scripts/general_scripts/20220504_MMEJ_detection_pipline.py "$PATH_TO_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/scripts/general_scripts/parsevcf.py "$PATH_TO_REALINGED_OUTPUT_DATA"/scripts
cp /home/labs/alevy/guyta/guy_master_project/scripts/general_scripts/pre_detection_processing.py "$PATH_TO_OUTPUT_DATA"/scripts
# cp 20220503_get_faste_context.sh "$PATH_TO_OUTPUT_DATA"/scripts
# cp /home/labs/alevy/guyta/guy_master_project/scripts/general_scripts/20220503_merge_VCF_with_reference_context.py "$PATH_TO_OUTPUT_DATA"/scripts

echo "This folder contains the data after processing with parsevcf.py" > ""$PATH_TO_REALINGED_OUTPUT_DATA"/docs/README.md"
python3 parsevcf.py -v ""$PATH_TO_DATA"/"$DATA_FILE_NAME"" -o ""$PATH_TO_REALINGED_OUTPUT_DATA"/"$DATA_FILE_NAME"_realined.vcf" -r $REF_GENOME_FILE

# Non-realigned
# python3 pre_detection_processing.py -v "$PATH_TO_DATA"/"$DATA_FILE_NAME" -w $CONTEXT_WINDOW_SIZE -r $REF_GENOME_FILE -o $PATH_TO_OUTPUT_DATA -d $FULL_DATE -sim "$SIMS" -ra "False"
# Realigned
python3 pre_detection_processing.py -v ""$PATH_TO_REALINGED_OUTPUT_DATA"/"$DATA_FILE_NAME"_realined" -w $CONTEXT_WINDOW_SIZE -r $REF_GENOME_FILE -o "$PATH_TO_REALINGED_OUTPUT_DATA" -d $FULL_DATE -sim "$SIMS" 

# Repair mechanisms (RM) detection
# Realigned
REALAIND_DATA_FOR_RM_DETECTION="$PATH_TO_REALINGED_OUTPUT_DATA"/"$FULL_DATE"_MMEJ_detection_formated_data.csv
python3 20220504_MMEJ_detection_pipline.py -i $REALAIND_DATA_FOR_RM_DETECTION -w $CONTEXT_WINDOW_SIZE -o "$PATH_TO_REALINGED_OUTPUT_DATA" -d $FULL_DATE -f $DATA_FILE_NAME -sim "$SIMS" 

# creating a README.md file for the pipline run
OUTPUT_DESCRIPTION="RM- repair mechanism\n# "$FULL_DATE"_"$CONTEXT_WINDOW_SIZE"_bp_context_window_"$DATA_FILE_NAME".bed -> A .bed file with coordinates of the indels +- context window\n# "$FULL_DATE"_getfasta_"$CONTEXT_WINDOW_SIZE"-> A folder with the reference genome context window\n# "$FULL_DATE"__MMEJ_detection_formated_data.csv -> A file with the data in the right format for the MMEJ detection process\n# "$FULL_DATE"_"$DATA_FILE_NAME"_classified_repair_mechanism.csv -> A dataframe with ALL the RM detection results and parameters\n# "$FULL_DATE"_"$DATA_FILE_NAME"_classified_repair_mechanism_prob.csv -> A dataframe with the RM detection results, only probabilities\n"
echo -e "This folder contains run of 20220503_main_pipeline.sh "$FULL_DATE" on the data from "$PATH_TO_DATA"\n INPUT data:\n"$DATA_FILE_NAME"\nOUTPUT:\n"$OUTPUT_DESCRIPTION"\n" >> "$PATH_TO_OUTPUT_DATA/docs/README.md"
echo "Done, output is at: "$PATH_TO_OUTPUT_DATA""


# bsub -q new-medium -R "rusage[mem=20000]" ./20220503_main_pipline.sh /home/labs/alevy/guyta/guy_master_project/results/general_pipline_results_test/20220505 /home/labs/alevy/guyta/guy_master_project/data/from_fabrizio/indel_sumlations/simulations test2


# bsub -q new-medium -R "rusage[mem=20000]" ./20220503_main_pipline.sh /home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20220512 /home/labs/alevy/guyta/guy_master_project/data/from_fabrizio/indel_sumlations/simulations test1 /home/labs/alevy/guyta/guy_master_project/data/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta 1745_20220512

# bsub -q new-medium -R "rusage[mem=20000]" ./20220503_main_pipline.sh /home/labs/alevy/guyta/guy_master_project/results/indels_simulations/simulations_from_Fabrizio/20220515 /home/labs/alevy/guyta/guy_master_project/data/from_fabrizio/indel_sumlations/simulations test1 /home/labs/alevy/guyta/guy_master_project/data/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta 1547_20220515

# TODO
# redocument the functions in ExpectationMaximization.py
# redecument the pre_processing script