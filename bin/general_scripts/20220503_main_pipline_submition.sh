# 2022.05.04
# defining flags to a bash script
while getopts 'h:f:v:o:m:' opt
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
    o)
        path_to_output="$OPTARG"
        ;;
    m)
        memory_required="$OPTARG"
        ;;
    esac
done

shift "$(( OPTIND - 1 ))"

# This script submits 20220503_main_pipline.sh as a job to the Wexac cluster
# module load samtools/1.9
# creating a README for the run
FULL_DATE=$(date +%H%M_%Y%m%d)
ARGS="$@"
echo -e "This folder contains the outputs of: $0\nLocated at: $(pwd) \nTime of submission: $FULL_DATE\nThe command was: $0 $ARGS" > "$path_to_output"/"$FULL_DATE"_README.txt
# indexing reference genome
# samtools faidx $reference_genome
# extracting file names
FILES=$(ls $path_to_vcfs/*.vcf | awk -F"/" '{print $NF}')
# iterating over the files and submiting each one of them as a seperate job
for file in $FILES;
do
    file_corrected=${file::-4}
    echo $file_corrected
    bsub -q new-short -m public_hosts -R "rusage[mem=$memory_required]" -J $file_corrected -e "$path_to_output"/"$FULL_DATE"_main_pipline_output_"$file_corrected"/err.txt -o "$path_to_output"/"$FULL_DATE"_main_pipline_output_"$file_corrected"/standart_output.txt ./20220503_main_pipline.sh $path_to_output $path_to_vcfs $file_corrected $reference_genome $FULL_DATE 1500 
    echo "Exit code $file: $?" >> "$path_to_output"/"$FULL_DATE"_README.txt 2>&1
done


# command example: 
# ./20220503_main_pipline_submition.sh -v /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/20220620_sims/mixed_files -o /home/labs/alevy/guyta/guy_master_project/results/indels_simulations/20220620_sims/mixed_files/20220703_3Kbp_res -f /home/labs/alevy/guyta/guy_master_project/data/indels_simulations/from_fabrizio/indel_sumlations/simulations/GWHAAEV00000000.1.genome.fasta -m 5000
