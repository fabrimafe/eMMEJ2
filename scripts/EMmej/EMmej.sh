#!/bin/bash

# defining flags to a bash script
# setting defult values
steps=2
mechanism_w=300
MM_w=1000

while getopts 'h:f:v:n:O:s:w:aci' opt
do 
    case "$opt" in
    h|\?)
        echo -e 'available options: \n-h : Argumants docs  \n-f : indexed (using tabix) reference genome (.fasta) \n-v : input vcf file (.vcf) \n-o : output path (string)\n-w : search window for pattern detection (int)\n-a: indication on whether an ancestral state is available\n-c: Indicate algorithm use on CRISPR data\n-s: indication of what steps to perform:\n  # 0: perform repair mechanism motif detection\n  # 1: perform repair mechanism motif detection & use Markov model\n  # 2: perform repair mechanism motif detection & use Markov model & estimate parameters using EM (int)\n-i: indicate wether to include context in output'
        exit 0
        ;;
    f)
        reference_genome="$OPTARG"
        ;;
    v)
        input_vcf="$OPTARG"
        ;;
    n)
        input_vcf_name="$OPTARG"
        ;;
    O)
        output_path="$OPTARG"
        ;;
    s)
        steps="$OPTARG"
        ;;
    w)
        mechanism_w="$OPTARG"
        ;;
    a)
        ancestral_state=-anc
        ;;
    c)
        crispr=-cr
        ;;
    i)
        includ_context=-ic
        ;;
    esac
done
ARGS="$@"
shift "$(( OPTIND - 1 ))"

FULL_DATE=$(date +%H%M_%Y%m%d)
output_dir="$output_path/$FULL_DATE"_"$input_vcf_name"_EMmej_run
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

echo -e "### Log ###\nStarting time: $(date +%H:%M-%d-%m-%Y)" >> "$output_dir"/README.txt
if [[ $steps == 0 ]] ;then
    # From raw vcf file to a dataframe with recognized repair mechanisms patterns
    python3 RMdetector.py -v "${input_vcf}/${input_vcf_name}.vcf" -o "${output_dir}/output_files" -r ${reference_genome} -w ${mechanism_w} ${crispr} ${ancestral_state} ${includ_context}
    echo "RMdetector DONE at: $(date +%H:%M-%d-%m-%Y)" ; echo "RMdetector DONE at: $(date +%H:%M-%d-%m-%Y)" >> "$output_dir"/README.txt
fi
if [[ $steps == 1 ]] ;then
    # From raw vcf file to a dataframe with recognized repair mechanisms patterns
    python3 RMdetector.py -v "${input_vcf}/${input_vcf_name}.vcf" -o "${output_dir}/output_files" -r ${reference_genome} -w ${mechanism_w} ${crispr} ${ancestral_state} ${includ_context}
    echo "RMdetector DONE at: $(date +%H:%M-%d-%m-%Y)" ; echo "RMdetector DONE at: $(date +%H:%M-%d-%m-%Y)" >> "$output_dir"/README.txt

    # From a dataframe with recognized repair mechanisms patterns to a dataframe with probability mesurments calculated using the Markov model
    python3 Markov_model_operator.py -v "$output_dir/output_files/RMdetector_output.tsv" -o "$output_dir/output_files" -w $MM_w -r $reference_genome ${includ_context}
    echo "Markov_model_operator DONE at: $(date +%H:%M-%d-%m-%Y)" ; echo "Markov_model_operator DONE at: $(date +%H:%M-%d-%m-%Y)" >> "$output_dir"/README.txt
fi
if [[ $steps == 2 ]] ;then
    # From raw vcf file to a dataframe with recognized repair mechanisms patterns
    python3 RMdetector.py -v "${input_vcf}/${input_vcf_name}.vcf" -o "${output_dir}/output_files" -r ${reference_genome} -w ${mechanism_w} ${crispr} ${ancestral_state} ${includ_context}
    echo "RMdetector DONE at: $(date +%H:%M-%d-%m-%Y)" ; echo "RMdetector DONE at: $(date +%H:%M-%d-%m-%Y)" >> "$output_dir"/README.txt
    # From a dataframe with recognized repair mechanisms patterns to a dataframe with probability mesurments calculated using the Markov model
    python3 Markov_model_operator.py -v "$output_dir/output_files/RMdetector_output.tsv" -o "$output_dir/output_files" -r $reference_genome -w $MM_w ${includ_context}
    echo "Markov_model_operator DONE at: $(date +%H:%M-%d-%m-%Y)" ; echo "Markov_model_operator DONE at: $(date +%H:%M-%d-%m-%Y)" >> "$output_dir"/README.txt
    # From a dataframe with probability mesurments calculated using the Markov model to an estimated indel length distributions and proportions using the EM
    python3 EM_operator.py -v "$output_dir/output_files/MM_output.tsv" -o "$output_dir/output_files" -r $reference_genome -w $mechanism_w # OPTIONAL ARGS: -b (int) -bs (float)
    echo "EM_operator DONE at: $(date +%H:%M-%d-%m-%Y)" ; echo "EM_operator DONE at: $(date +%H:%M-%d-%m-%Y)" >> "$output_dir"/README.txt
    
fi

END_TIME=$(date +%H:%M-%d-%m-%Y)
echo -e "### EMmej DONE ###\nTime: $END_TIME" >> "$output_dir"/README.txt
echo -e "### EMmej DONE ###\nTime: $END_TIME \nOutput is at:\n$output_dir"

conda deactivate

# command example:
# ./EMmej.sh -v /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/20220620_sims/mixed_files/20220629_07_MMEJ_03_NHEJ_mix_unlabled.vcf -o /home/labs/alevy/guyta/guy_master_project/results/EMmej_tests -f /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta -s 2 -c True
# ./EMmej_submition_to_cluster.sh -v /home/labs/alevy/guyta/guy_master_project/results/soybean/Liu_et.al.2020/annotation_intersection/20220723_chr01_04_06_08_09_EMmej_format_1000_indels_samples/20220723_chr01_04_06_08_09_EMmej_exon -o /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/Mutation_bias_s41586-021-04269-6/EMmej/20220723_exon_only -f /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta -m 3000

# example with crispr data:
# bsub -q new-short -m public_hosts -R "rusage[mem=3000]" -J NHEJ_and_MMEJ -e /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/Mutation_bias_s41586-021-04269-6/EMmej/20220814_whole_genome/err.txt -o /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/Mutation_bias_s41586-021-04269-6/EMmej/20220814_whole_genome/stdout.txt ./EMmej.sh -v /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/Mutation_bias_s41586-021-04269-6/raw_data_analyzing/raw_data_from_Grey_res -n 20220717_high_confidance_set_indels_EMmej_format -O /home/labs/alevy/guyta/guy_master_project/results/arabidopsis/Mutation_bias_s41586-021-04269-6/EMmej/20220814_whole_genome -f /home/labs/alevy/guyta/guy_master_project/data/arabidopsis/TAIR10_ref_genome/TAIR10_fasta/GCF_000001735.4_TAIR10.1_genomic.fna -a False -s 2 -c True

# example on test data
# ./EMmej.sh -v /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/test_data -n EMmej_test_data_indels -O /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/test_results -f /home/labs/alevy/guyta/guy_master_project/tests/EMmej_tests/test_data_reference_genome/EMmej_test_data_ref_genome.fa -s 0 -a -c -i

