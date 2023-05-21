#!/bin/bash

module load bedtools

VCF='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/raw_data_to_EMmej_format/20221126/DGRP_EMmej_formated_data.vcf.gz'
OUT_PTH='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/DGRP/annotation_intersection/20221126'
REP_ANN='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/dm6_reptet_ann_chr_renamed_sorted_merged.bed'
GENE_ANN='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/GCF_000001215.4_genic_merged_sorted.bed.gz'
EXON_ANN='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/GCF_000001215.4_exons_merged_sorted.bed.gz'

# Repeat and non-repeat
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${VCF} )) -b ${REP_ANN} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_repeated_regions.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${VCF} )) -b ${REP_ANN} -wa -f 1.0 -header -v | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions.vcf.gz"


# Exons and non genic
NON_REP_DATA="${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b ${EXON_ANN} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_exons.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b ${GENE_ANN} -wa -f 1.0 -header -v | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_non_genic.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b ${GENE_ANN} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_genic.vcf.gz"


# Histon modifications:
HIST_MOD_ANN_PTH='/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/table4026.gff3_Adult_femail_hist_mods_stratified/'
# AdultFemale_8wg16_ChIP_chip
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_8wg16_ChIP_chip_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_8wg16_ChIP_chip.vcf.gz"
# AdultFemale_CBP_ChIP_chip
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_CBP_ChIP_chip_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_CBP_ChIP_chip.vcf.gz"
# AdultFemale_CBP_chip_seq_merged
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_CBP_chip_seq_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_CBP_chip_seq.vcf.gz"
# H3K27Ac_chip_seq
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_H3K27Ac_chip_seq_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_H3K27Ac_chip_seq.vcf.gz"
# H3K4Me1_chip_seq
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_H3K4Me1_chip_seq_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_H3K4Me1_chip_seq.vcf.gz"
# H3K4Me3_chip_seq
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_H3K4Me3_chip_seq_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_H3K4Me3_chip_seq.vcf.gz"
# H3K9Ac_chip_seq
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_H3K9Ac_chip_seq_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_H3K9Ac_chip_seq.vcf.gz"
# H3K9AC
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_ID=AdultFemale_H3K9AC_merged.gff3" -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_H3K9AC.vcf.gz"
# non histon modifications
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${NON_REP_DATA} )) -b "${HIST_MOD_ANN_PTH}/table4026.gff3_AdultFemale_merged.gff3" -v -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_non_repeated_regions_AdultFemale_no_histon_modifications.vcf.gz"

# repeated regions stratified to famalies:
REP_DATA="${OUT_PTH}/DGRP_EMmej_formated_repeated_regions.vcf.gz"

DNA_TE_ANNOT="/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/dm6_repeat_ann_stratification/dm6_DNA_TE_sorted_merged.bed"
LINE_TE_ANNOT="/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/dm6_repeat_ann_stratification/dm6_LINE_TE_sorted_merged.bed"
LTR_TE_ANNOT="/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/dm6_repeat_ann_stratification/dm6_LTR_TE_sorted_merged.bed"
RC_TE_ANNOT="/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/dm6_repeat_ann_stratification/dm6_RC_TE_sorted_merged.bed"
SATELLITE_ANNOT="/home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/annotations/dm6_repeat_ann_stratification/dm6_Satellite_sorted_merged.bed"

bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${REP_DATA} )) -b ${DNA_TE_ANNOT} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_DNA_TE.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${REP_DATA} )) -b ${LINE_TE_ANNOT} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_LINE_TE.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${REP_DATA} )) -b ${LTR_TE_ANNOT} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_LTR_TE.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${REP_DATA} )) -b ${RC_TE_ANNOT} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_RC_TE.vcf.gz"
bedtools intersect -a <(awk '{print $1"\t"$2"\t"$2+length($3)"\t"$3"\t"$4"\t"$5}' <( zcat ${REP_DATA} )) -b ${SATELLITE_ANNOT} -wa -f 1.0 -header | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gzip > "${OUT_PTH}/DGRP_EMmej_formated_Satelite.vcf.gz"