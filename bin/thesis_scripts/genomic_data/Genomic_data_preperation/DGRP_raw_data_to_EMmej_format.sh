#!/bin/bash

module load bcftools
module load vcftools

RAW_DATA='/home/labs/alevy/guyta/guy_master_project/data/Drosophila/DGRP/vcfs/dgrp2.vcf'
OUT_PTH='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/raw_data_to_EMmej_format/20221126'
MIN_MISS=0.02
MAC=2
# filtering for indels only
bcftools view --types indels -Ov \
    "${RAW_DATA}"  | gzip > "${OUT_PTH}/DGRP_Indels.vcf.gz"  # | head -n 10000


# Transforming multi to mono allelic format
cat "${OUT_PTH}/DGRP_Indels.vcf.gz" | bcftools norm -Ov -m-any | gzip > "${OUT_PTH}/norm_DGRP_Indels.vcf.gz"

# filtering for indels with at least ${MIN_MISS} data availability
vcftools --gzvcf "${OUT_PTH}/norm_DGRP_Indels.vcf.gz" --max-missing ${MIN_MISS} --recode --recode-INFO-all --out "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}" 
cat "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}.recode.vcf" | gzip > "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}.recode.vcf.gz"

# filtering for indels with at least ${MAC} minor allele count
vcftools --gzvcf "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}.recode.vcf.gz" --mac ${MAC} --recode --recode-INFO-all --out "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}_MAC${MAC}" 
cat "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}_MAC${MAC}.recode.vcf" | gzip > "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}_MAC${MAC}.recode.vcf.gz"

# calculating allele freq
vcftools --gzvcf "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}_MAC${MAC}.recode.vcf.gz" --freq --out "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}_MAC${MAC}_AF" 

module load miniconda
conda activate guy_mmej_env
# Using the allele freq to set ancestral stat
python3 add_ancesral_state_column.py -v "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}_MAC${MAC}.recode.vcf.gz" -af "${OUT_PTH}/norm_DGRP_Indels_min_missing_data_${MIN_MISS}_MAC${MAC}_AF.frq" -o ${OUT_PTH}

conda deactivate


cat "${OUT_PTH}/DGRP_EMmej_formated_data.vcf" |  gzip > "${OUT_PTH}/DGRP_EMmej_formated_data.vcf.gz"

rm "${OUT_PTH}/DGRP_EMmej_formated_data.vcf"

echo -e "This folder contains the peocess of cleaning the data from: \n - SNPs\n - Indels that present in less the ${MIN_MISS} of the accessions\n - Indels with minor allele count < ${MAC}\n Than we calculate allele frequency and generate the ancestral state column and the EMmej format.\n script is: \n/home/labs/alevy/guyta/guy_master_project/scripts/Drosophila/Clark_AG_et_al_G3_2015/construct_ancestral_state/Clark_raw_data_to_EMmej_format.sh" > "${OUT_PTH}/README.md"


