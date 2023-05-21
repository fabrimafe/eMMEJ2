#!/bin/bash

# This script submmit a EMmej job to the wexac cluster

# 2022.05.04
# defining flags to a bash script
while getopts 'h:f:v:O:m:ca' opt
do 
    case "$opt" in
    h|\?)
        echo -e 'available options: -h\n-f : reference genome in a fasta format\n-v : a path to the vcf files\n-o : output path\n-m : memory required (RAM in Mb)'
        ;;
    f)
        reference_genome="$OPTARG"
        ;;
    v)
        path_to_vcfs="$OPTARG"
        ;;
    O)
        path_to_output="$OPTARG"
        ;;
    m)
        memory_required="$OPTARG"
        ;;
    c)
        crispr=-cr
        ;;
    a)
        ancestral_state=-a
        ;;
    esac
done

shift "$(( OPTIND - 1 ))"

# This script submits 20220503_main_pipline.sh as a job to the Wexac cluster
module load samtools/1.9
# creating a README for the run
FULL_DATE=$(date +%H%M_%Y%m%d)
ARGS="$@"
# echo -e "This folder contains the outputs of: $0\nLocated at: $(pwd) \nTime of submission: $FULL_DATE\nThe command was: $0 $ARGS" > "$path_to_output"/"$FULL_DATE"_README.txt
# indexing reference genome
samtools faidx $reference_genome
# extracting file names
FILES=$(ls $path_to_vcfs/*.vcf | awk -F"/" '{print $NF}')
# iterating over the files and submiting each one of them as a seperate job
module load miniconda
conda activate guy_mmej_env
for file in $FILES;
do
    file_corrected=${file::-4}
    echo $file_corrected
    # echo "$path_to_vcfs/$file_corrected"
    bsub -q new-short -e /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_err_%J.txt -o /home/labs/alevy/guyta/guy_master_project/results/my_err_and_output/std_output_%J.txt -m public_hosts -R "rusage[mem=$memory_required]" -J $file_corrected  "python3 EMmej.py -v $path_to_vcfs/$file_corrected.vcf -o $path_to_output -r $reference_genome $ancestral_state $crispr -vr -p markov -db -e 12 -wm 150" #  -wm 50 -wp 50 -cr -wm 75 -b 100   

    # echo "Exit code $file: $?" >> "$path_to_output"/"$FULL_DATE"_README.txt 2>&1
    # sleep 3s
    # mv ""${path_to_output}"/"${FULL_DATE}"_README.txt" ""${path_to_output}/${FULL_DATE}"_"${file_corrected}"_EMmej_run/"
done


# command example: 
# ./EMmej_submition_to_cluster.sh -f /home/labs/alevy/guyta/guy_master_project/data/arabidopsis/TAIR10_ref_genome/TAIR10_fasta/GCF_000001735.4_TAIR10.1_genomic.fna -v /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/Mutation_bias_s41586-021-04269-6/raw_data_analyzing/raw_data_from_Grey_res -O /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/Mutation_bias_s41586-021-04269-6/EMmej/20220903_whole_genome -m 3000 -c True
# ./EMmej_submition_to_cluster.sh -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/stratification_by_annotation/block_sampaled_vcfs/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/20221010_5000000bp_blocks_bootN_5 -O /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/20221010 -f /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa -m 2000 -a
# ./EMmej_submition_to_cluster.sh -v /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/stratification_by_annotation/block_sampaled_vcfs/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/20221016_2000000bp_blocks -O /home/labs/alevy/guyta/guy_master_project/results/Drosophila/Clark_AG_et_al_G3_2015/EMmej_results/GDL_indel_EMmej_format_non_repeated_regions_exon_no_MT/20221021 -f /home/labs/alevy/guyta/guy_master_project/results/Drosophila/ref_genomes/D_melanogaster/ref_genome/ISO1_GCF_000001215.4.fa -m 2000 -a