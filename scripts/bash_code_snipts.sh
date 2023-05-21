#!/bin/bash
ARGS="$@" # get all arguments that was passed after the script's name
# an example for a README format (there is an additional Exit code logging in the last command of this script):
FULL_DATE=$(date +%H%M_%Y%m%d)
echo -e "This folder contains the outputs of: $0\nLocated at: $(pwd) \nTime of submission: $FULL_DATE\nThe command was: $0 $ARGS" > logfile.txt


echo $0 # name of the script
echo "$#" # number of virables that were passed to the sctipt
# example to a for loop with if statments
vars=$( echo 1{4,5}{0..9} )
for v in $vars
do
    if [[ $v > 150 ]]
    then
        echo "$v is greater then 150"
    elif [[ $v = 150 ]]
    then
        echo "$v is equal to 150"
    elif [[ $v < 150 ]]
    then
        echo "$v is less then 150"
    fi
done

vars=$( echo 0{1,9} )
for v in $vars
do
    echo $v
done



# defining flags to a bash script
while getopts 'h:f:v:o:' opt
do 
    case "$opt" in
    h|\?)
        echo 'available options: -h  -f [filename] -v [vcffile -o [outputpath]'
        ;;
    f)
        input_fasta="$OPTARG"
        ;;
    v)
        input_vcf="$OPTARG"
        ;;
    o)
        output_path="$OPTARG"
        ;;
    esac
done

shift "$(( OPTIND - 1 ))"

echo $input_fasta
echo $input_vcf
echo $output_path

# read arg # ask for input
# echo $arg > stout.txt # saves input into a file
# read next_arg # ask for another input
# echo $next_arg >> stout.txt # saves input into a file without overiting

# adding the Exit code to logfile.txt
echo "Exit code: $?" >> logfile.txt 2>&1

