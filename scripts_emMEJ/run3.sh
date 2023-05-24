#!/bin/bash
module load miniconda
module load tabix
module load bedtools
conda activate emmej

reference=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_arabidopsis/data/Arabidopsis_thaliana.fa #/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/data/Drosophila_melanogaster.fa
xfolder=/home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Arabidopsis1001/filtered/mini #/home/labs/alevy/fabrizio/workspace/guy/rerun/03302023/Clark/subsampled1k/16/mini
echo "run RM"
python ./scripts/EMmej/RMdetector.py --vcf ${xfolder}/indels.vcf --ref $reference --anc -o ${xfolder}/RMoutput.tsv -w 50 -ic #-mhl 2
echo "MM"
#python ./scripts/EMmej/Markov_model_operator.py --vcf ${xfolder}/RMoutput.tsv --ref $reference -p markov -o ${xfolder}/MM_output.tsv

