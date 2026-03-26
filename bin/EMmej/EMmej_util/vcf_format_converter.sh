#!/bin/bash

# This script convert any vcf format to the right format for EMmej

# defining flags to a bash script
while getopts 'h:r:v:f:O:' opt
do 
    case "$opt" in
    h|\?)
        echo -e 'available options: \n-h : Argumants docs  \n-r : indexed (using tabix) reference genome (.fa & .fai) \n-v : input vcf file (.vcf) \n-o : output path (string)\n-f: indication whether an AF file (vcftools --freq2 .frq format) is available (bool)'
        exit 0
        ;;
    r)
        reference_genome="$OPTARG"
        ;;
    v)
        input_vcf="$OPTARG"
        ;;
    O)
        output_path="$OPTARG"
        ;;
    f)
        allele_freq="$OPTARG"
        ;;
    esac
done
ARGS="$@"
shift "$(( OPTIND - 1 ))"

module load vcftools

vcftools --vcf ${input_vcf} --freq2 --out "${output_path}/allele_freq"

module load miniconda
conda activate guy_mmej_env

# python3 vcf_format_converter.py -v ${input_vcf} -o ${output_path} -r ${reference_genome} -af ${allele_freq}
python3 vcf_format_converter.py -v ${input_vcf} -o ${output_path} -r ${reference_genome} -af "${output_path}/allele_freq.frq"

conda deactivate



# Command example:
# ./vcf_format_converter.sh -r /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa  -v /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/GDL_Indels_reformated.vcf -O /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter -f /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/vcf_converter/allele_freq.frq
