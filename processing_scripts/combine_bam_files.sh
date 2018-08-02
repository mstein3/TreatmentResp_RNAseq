#!/bin/bash

# Script to combine .bam files with the same Extract Order ID into a single .bam file
# Will let verifyBAM work better. 
# the read command is going to look at whatever file you pipe into the script with the < operator

INPUT_DIR=/group/nicolae-lab/users/mstein3/STARhg19_May2017_bamFiles
OUTPUT_DIR=/group/nicolae-lab/users/mstein3/STARhg19_CombinedBAM

while read LINE; do
   	OUTPUT_NAME="$(basename $LINE)_STARhg19_combined.bam"
   	cat $INPUT_DIR/$LINE* > $OUTPUT_DIR/$OUTPUT_NAME
    echo "Combining $LINE .bam files and write to /group/nicolae-lab/users/mstein3/STARhg19_CombinedBAM/$NEWFILE"
done < $1

#done < /group/ober-resources/users/mstein3/rna.seq.2017/processing_scripts/uniq_numbers_of_bam_files.txt


# This .sh would take a long time, might be safer to run as a .pbs script. 
