#!/bin/bash
module load miniconda
module load tabix
module load bedtools
conda activate emmej
for xdataset in DGRP Clark Arabidopsis1001;do
outroot=/home/labs/alevy/fabrizio/workspace/guy/rerun/03302023
if [ $xdataset == "Arabidopsis1001" ];then
reference=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_arabidopsis/data/Arabidopsis_thaliana.fa
filter=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_arabidopsis/results/toast/referenceancestral_map35.bed.gz
inputvcf=/home/labs/alevy/guyta/guy_master_project/results/arabidopsis/1001_genomes/raw_data_analyzing/raw_data_to_EMmej_format_mac2_no_1bp_del/1001genomes_EMmej_formated_data.vcf.gz
elif [ $xdataset == "DGRP" ];then
reference=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/data/Drosophila_melanogaster.fa
filter=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/results/toast/referenceancestral_map35.bed.gz
inputvcf=/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/raw_data_to_EMmej_format/20221126/DGRP_EMmej_formated_data.vcf.gz
elif [ $xdataset == "Clark" ];then
reference=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/data/Drosophila_melanogaster.fa
filter=/home/labs/alevy/fabrizio/workspace/ancestralstate/msa_dro/results/toast/referenceancestral_map35.bed.gz
inputvcf=/home/labs/alevy/fabrizio/workspace/guy/rerun/Drosophila_clark.vcf.gz
fi
echo $xdataset
xfolder=${outroot}/${xdataset}
mkdir -p $xfolder
zcat $inputvcf | head -1 > ${xfolder}/indels.vcf
#bedtools intersect -b $filter -a <( zcat $inputvcf | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5}' ) | awk -v OFS='\t' '{print $1,$3,$4,$5,$6}' >> ${xfolder}/indels.vcf
cd /home/labs/alevy/fabrizio/github/guymasterproject_final/guy2023
echo "RM"
python ./scripts/EMmej/RMdetector.py --vcf ${xfolder}/indels.vcf --ref $reference --anc -o ${xfolder}/RMoutput.tsv
echo "MM"
python ./scripts/EMmej/Markov_model_operator.py --vcf ${xfolder}/RMoutput.tsv --ref $reference -p markov -o ${xfolder}/MM_output.tsv
echo "EM"
python ./scripts/EMmej/EM_operator.py --vcf ${xfolder}/MM_output.tsv --ref $reference -o ${xfolder}/EM_output.tsv
done

