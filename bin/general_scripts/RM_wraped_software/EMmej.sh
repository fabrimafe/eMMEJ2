#!/bin/bash

# defining flags to a bash script
while getopts 'h:f:v:o:s:' opt
do 
    case "$opt" in
    h|\?)
        echo -e 'available options: \n-h : Argumants docs  \n-f : indexed (using tabix) reference genome (.fasta) \n-v : input vcf file (.vcf) \n-o : output path (string)\n-s: indication of what steps to perform:\n #0: perform repair mechanism motif detection\n #1: perform repair mechanism motif detection & use Markov model\n #2: perform repair mechanism motif detection & use Markov model & estimate parameters using EM (int)'
        exit 0
        ;;
    f)
        reference_genome="$OPTARG"
        ;;
    v)
        input_vcf="$OPTARG"
        ;;
    o)
        output_path="$OPTARG"
        ;;
    s)
        steps="$OPTARG"
        ;;
    esac
done

shift "$(( OPTIND - 1 ))"

ARGS="$@"
FULL_DATE=$(date +%H%M_%Y%m%d)
output_dir="$output_path/$FULL_DATE"_EMmej_run
# Backup run's scripts
mkdir $output_dir
mkdir $output_dir/scripts
mkdir $output_dir/docs
mkdir $output_dir/output_files
scripts_dir=$( pwd )
modules_dir="/home/labs/alevy/guyta/guy_master_project/src/"
cp "$scripts_dir"/EMmej.sh "$scripts_dir"/RMdetector.py "$scripts_dir"/Markov_model_operator.py "$scripts_dir"/EM_operator.py "$output_dir"/scripts
cp "$modules_dir"/MicroHomology_module_v3.py "$modules_dir"/MMEJ_2nd_order_MM_v2.py "$modules_dir"/ExpectationMaximization_q.py "$modules_dir"/pysam_getfasta.py "$output_dir"/scripts

echo -e "This folder contains the outputs of: $0\nLocated at: $(pwd) \nTime of submission: $FULL_DATE\nThe command was: $0 $ARGS\n\nAll scripts that got involved in the process are at: scripts\nAll output files are at: output_file\nAll logs and documentations are at: docs" > "$output_dir"/README.txt

module load miniconda
conda activate guy_mmej_env

if [[ $steps == 0 ]] ;then
    # From raw vcf file to a dataframe with recognized repair mechanisms patterns
    python3 RMdetector.py -v $input_vcf -o "$output_dir/output_files" -w 300 -ra True -anc False -r $reference_genome
    echo "RMdetector.py DONE"
fi
if [[ $steps == 1 ]] ;then
    # From raw vcf file to a dataframe with recognized repair mechanisms patterns
    python3 RMdetector.py -v $input_vcf -o "$output_dir/output_files" -w 300 -ra True -anc False -r $reference_genome
    echo "RMdetector.py DONE"
    # From a dataframe with recognized repair mechanisms patterns to a dataframe with probability mesurments calculated using the Markov model
    python3 Markov_model_operator.py -v "$output_dir/output_files/RMdetector_output.csv" -o "$output_dir/output_files" -w 1500 -r $reference_genome -ir False
    echo "Markov_model_operator.py DONE"
fi
if [[ $steps == 2 ]] ;then
    # From raw vcf file to a dataframe with recognized repair mechanisms patterns
    python3 RMdetector.py -v $input_vcf -o "$output_dir/output_files" -w 300 -ra True -anc False -r $reference_genome
    echo "RMdetector.py DONE"
    # From a dataframe with recognized repair mechanisms patterns to a dataframe with probability mesurments calculated using the Markov model
    python3 Markov_model_operator.py -v "$output_dir/output_files/RMdetector_output.csv" -o "$output_dir/output_files" -w 1500 -r $reference_genome -ir False
    echo "Markov_model_operator.py DONE"
    # From a dataframe with probability mesurments calculated using the Markov model to an estimated indel length distributions and proportions using the EM
    python3 EM_operator.py -v "$output_dir/output_files/MM_output.csv" -o "$output_dir/output_files" -w 15 -r $reference_genome -b 10 -bs 0.7
    echo "EM_operator.py DONE"
fi

# command example:
# ./EMmej.sh -v /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/20220620_sims/mixed_files/20220629_07_MMEJ_03_NHEJ_mix_unlabled.vcf -o /home/labs/alevy/guyta/guy_master_project/scripts/general_scripts/RM_wraped_software/EMmej_test -f /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta