#!/bin/bash
module load miniconda
module load tabix
module load bedtools
conda activate emmej

#cat ~/workspace/guy/simulations/Athaliana/original_sims_EMmej_format/deletion_MMEJ_l2_d50_EMmej_format.vcf | head -101 > ~/workspace/guy/simulations/Athaliana/original_sims_EMmej_format/mini/indels.vcf
reference=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_arabidopsis/data/Arabidopsis_thaliana.fa #/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/data/Drosophila_melanogaster.fa
xfolder=~/workspace/guy/simulations/Athaliana/original_sims_EMmej_format/mini
#xfolder=/home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Arabidopsis1001/filtered/mini #/home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Clark/subsampled1k/16/mini
echo "run RM"
python ./scripts/EMmej/RMdetector.py --vcf ${xfolder}/indels.vcf --ref $reference -o ${xfolder}/RMoutput.tsv -w 50 -ic -mhl 4,3,2,1 #3,2,1 #2,1,0 #--anc
echo "MM"
#python ./scripts/EMmej/Markov_model_operator.py --vcf ${xfolder}/RMoutput.tsv --ref $reference -p markov -o ${xfolder}/MM_output.tsv
echo "EM"
#python ./scripts/EMmej/EM_operator.py --vcf ${xfolder}/RMoutput.tsv --ref $reference -o ${xfolder}/EM_output.tsv -e 1 -mhl 4 -ildt "savitzky_golay" #"uniform" -mhl 2,3,4,0
#python ./scripts/EMmej/EM_operator.py --vcf ${xfolder}/RMoutput.tsv --ref $reference -o ${xfolder}/EM_output.tsv -e 1 -mhl 0,2 -ildt "savitzky_golay" -dpol #-c 0.000000001
python ./scripts/EMmej/EM_operator.py --vcf ${xfolder}/RMoutput.tsv --ref $reference -o ${xfolder}/EM_output.tsv -e 1 -mhl 4,3,2,1 -ildt "savitzky_golay" #-ildt "full" # -dpol #-c 0.000000001



