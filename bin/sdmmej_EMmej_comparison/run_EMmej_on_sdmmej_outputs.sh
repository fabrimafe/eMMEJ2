#!/bin/bash
# This script takes the EMmej VCF format + reference genomes produced by: /home/labs/alevy/guyta/guy_master_project/scripts/sdmmej_EMmej_comparison/run_EMmej_on_sdmmej_outputs.sh
# and run EMmej
cd /home/labs/alevy/guyta/guy_master_project/scripts/EMmej
FULL_DATE=$(date +%Y%m%d)
OUTPUT_PATH="/home/labs/alevy/guyta/guy_master_project/results/sdmmej_EMmej_comparison/${FULL_DATE}"

# mkdir "${OUTPUT_PATH}"
# # R0_WT
# ./EMmej_submition_to_cluster.sh -v "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_WT/20221008" -f "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_WT/20221008/20221008_R0_WT_sdmmej_output_EMmej_input_ref.fa" -O "${OUTPUT_PATH}" -m 1000 -c True -a False

# # R0_POLQ
# ./EMmej_submition_to_cluster.sh -v "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_POLQ/20221008" -f "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_POLQ/20221008/20221008_R0_POLQ_sdmmej_output_EMmej_input_ref.fa" -O "${OUTPUT_PATH}" -m 1000 -c True -a False

# # R0_Lig4
# ./EMmej_submition_to_cluster.sh -v "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_Lig4/20221008" -f "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_Lig4/20221008/20221008_R0_Lig4_sdmmej_output_EMmej_input_ref.fa" -O "${OUTPUT_PATH}" -m 1000 -c True -a False


mkdir "${OUTPUT_PATH}"
# R0_WT
./EMmej_submition_to_cluster.sh -v "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_WT/20221030" -f "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_WT/20221030/20221030_R0_WT_sdmmej_output_EMmej_input_ref.fa" -O "${OUTPUT_PATH}" -m 1000 -c True -a False

# R0_POLQ
./EMmej_submition_to_cluster.sh -v "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_POLQ/20221030" -f "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_POLQ/20221030/20221030_R0_POLQ_sdmmej_output_EMmej_input_ref.fa" -O "${OUTPUT_PATH}" -m 1000 -c True -a False

# R0_Lig4
./EMmej_submition_to_cluster.sh -v "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_Lig4/20221030" -f "/home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_Lig4/20221030/20221030_R0_Lig4_sdmmej_output_EMmej_input_ref.fa" -O "${OUTPUT_PATH}" -m 1000 -c True -a False


# mkdir "${OUTPUT_PATH}"
# module load miniconda
# conda activate guy_mmej_env
# # R0_WT
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=900]" "python3 EMmej.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_WT/20221030 -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_WT/20221030/20221030_R0_WT_sdmmej_output_EMmej_input_ref.fa -o ${OUTPUT_PATH} -p markov -db -cr -wm 75 -wp 75 -vr -e 12"

# # R0_POLQ
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=900]" "python3 EMmej.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_POLQ/20221030 -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_POLQ/20221030/20221030_R0_POLQ_sdmmej_output_EMmej_input_ref.fa -o ${OUTPUT_PATH} -p markov -db -cr -wm 75 -wp 75 -vr -e 12"

# # R0_Lig4
# bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=900]" "python3 EMmej.py -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_Lig4/20221030 -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Terrence_Hanscom_NAR_2022/sdmmej_EMmej_comparison/sdmmej_output_to_EMmej_input/R0_Lig4/20221030/20221030_R0_Lig4_sdmmej_output_EMmej_input_ref.fa -o ${OUTPUT_PATH} -p markov -db -cr -wm 75 -wp 75 -vr -e 12"
