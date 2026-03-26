#!/bin/bash

# This script takes a VCF file and produce a .bed file with the following format:
# CHR   Start   End
# based on indel coordinates and a given window size

PATH_TO_DATA=$1
DATA_FILE_NAME=$2
CONTEXT_WINDOW_SIZE=$3
PATH_TO_OUTPUT_DATA=$4
FULL_DATE=$5

cp 20220503_indels_context_window_coordination_generation.sh "$PATH_TO_OUTPUT_DATA"

awk -F "\t" -v window=$CONTEXT_WINDOW_SIZE '{ if($1 !~ "^#") {print $1"\t"$2-window"\t"$2+window}}' "$PATH_TO_DATA"/"$DATA_FILE_NAME".vcf > "$PATH_TO_OUTPUT_DATA"/"$FULL_DATE"_"$CONTEXT_WINDOW_SIZE"bp_context_window_"$DATA_FILE_NAME".bed

echo ""$PATH_TO_OUTPUT_DATA"/"$FULL_DATE"_"$CONTEXT_WINDOW_SIZE"bp_context_window_"$DATA_FILE_NAME".bed is the output of: 20220503_indels_context_window_coordination_generation.sh, $FULL_DATE" > "$PATH_TO_OUTPUT_DATA"/indels_context_README.md