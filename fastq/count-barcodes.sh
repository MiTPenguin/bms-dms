#!/usr/bin/env bash

# Usage: count-barcodes.sh <FASTQ_FILE> <OUTPUT_FILE> 

export FASTQ_FILE=$1
export OUTPUT_FILE=$2

zcat $FASTQ_FILE | \
    awk '{if(NR%4==2) print substr($0,1,21)}' | \
    sort | \
    uniq -c | \
    sed 's/^ *//g' | \
    sed 's/ /\t/g' > $OUTPUT_FILE
