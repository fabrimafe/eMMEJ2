#!/bin/bash


# 03.05.2022
# This script generates a .bed file containing fasta context that corespond to 
# coordinates from a file with the following format:
# Chr    start_pos   end_pos
# and a reference genome

# module load bedtools
DATE=$1
FULL_DATE=$2
CONTEXT_WINDOW_SIZE=$3
PATH_TO_DATA=$4
DATA_FILE_NAME=$5
REF_GENOME_FILE=$6
PATH_TO_OUTPUT_DATA=$7
PATH_TO_OUTPUT_DATA="$PATH_TO_OUTPUT_DATA"/"$FULL_DATE"_getfasta_$CONTEXT_WINDOW_SIZE
mkdir $PATH_TO_OUTPUT_DATA

echo "### getfasta started"
bedtools getfasta -fi $REF_GENOME_FILE -bed "$PATH_TO_DATA$DATA_FILE_NAME" > "$PATH_TO_OUTPUT_DATA"/"$DATA_FILE_NAME"

echo "This folder contains run of 20220503_get_faste_context.sh "$FULL_DATE" on the data from "$PATH_TO_DATA", INPUT data: "$DATA_FILE_NAME", OUTPUT: "$DATE"_"$DATA_FILE_NAME"_"$CONTEXT_WINDOW_SIZE"bp.bed -> A file containing all the fasta contexts that coresponding to the coordinates from the input data" > "$PATH_TO_OUTPUT_DATA"/README.md

echo "Done, output is at: "$PATH_TO_OUTPUT_DATA"/"$FULL_DATE"_getfasta_$CONTEXT_WINDOW_SIZE"
