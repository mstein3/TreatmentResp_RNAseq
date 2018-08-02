#!/bin/bash

# the read command is going to look at whatever file you pipe into the script with the < operator
while read LINE; do
    NEWFILE=$(echo "$LINE" | cut -f1)
    OLDFILEPATH=$(echo "$LINE" | awk '{print $7}')
    cp $OLDFILEPATH /group/ober-resources/users/mstein3/RNAseq/Reseq.All.fastq/$NEWFILE
    echo "Moving $OLDFILEPATH to /group/ober-resources/users/mstein3/RNAseq/Reseq.All.fastq/$NEWFILE"
done < /group/ober-resources/users/mstein3/RNAseq/FC.reseq.raw/RFSS_4_noheader.txt

