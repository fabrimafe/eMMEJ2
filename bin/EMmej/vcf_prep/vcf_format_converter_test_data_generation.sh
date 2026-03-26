#!/bin/bash
# 27.09.2022

# This folder contains all steps that I took in order to develop: vcf_converter.sh

# The 1st step is to generate some fake example dataset

module load bcftools

# Creating a sample dataset based on /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Clark_AG_et_al_G3_2015
bcftools view /home/labs/alevy/guyta/guy_master_project/data/Drosophila/Clark_AG_et_al_G3_2015/GDL_Indels.vcf | head -n 300 > /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_regular_format.vcf

# Converting to monoalelic format
bcftools norm -Ov -m-any /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_regular_format.vcf > /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_regular_monoalelic_format.vcf

awk -F "\t" '{ if($1!~ /^#/){ print }}' /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_regular_monoalelic_format.vcf > /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_regular_monoalelic_format_no_headers.vcf

module load miniconda
conda activate guy_mmej_env

python3 vcf_allel_format_converter_for_test_data_generation.py

conda deactivate

cat /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_regular_format.vcf | head -n 29 > /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_reformated.vcf

cat /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_reformated_no_headers.vcf >> /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_reformated.vcf