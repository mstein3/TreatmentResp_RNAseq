#!/bin/bash -x

# the read command is going to look at whatever file you pipe into the script with the < operator
while read LINE; do
		
        echo "qsub -v FILE=$LINE MakeBai_IndivBam.pbs"
        qsub -v FILE=$LINE MakeBai_IndivBam.pbs
        
done < $1


# add an extract order file as $1: in/group/ober-resources/users/mstein3/rna.seq.2017/processing_scripts/ExtractOrd_plain.txt
