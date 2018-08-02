#!/bin/bash -x

# the read command is going to look at whatever file you pipe into the script with the < operator
while read COL_ONE_BAM COL_TWO_FINDIV; do
		#echo "qsub ONE=$COL_ONE_BAM MakeBai_CombBam.pbs"
		#qsub -v ONE=$COL_ONE_BAM MakeBai_CombBam.pbs
        echo "qsub -v ONE=$COL_ONE_BAM,TWO=$COL_TWO_FINDIV verifyBAM.pbs"
        qsub -v ONE=$COL_ONE_BAM,TWO=$COL_TWO_FINDIV verifyBAM.pbs
        
done < $1

# add the smID file as $1. 
# smid information for the samtools MakeBai_CombBam.pbs: /group/ober-resources/users/mstein3/rna.seq.2017/processing_scripts/smID_for_samtools.txt
#smID information for verifyBamID: /group/ober-resources/users/mstein3/rna.seq.2017/processing_scripts/smID_ImHTgx_SortedVerifyBAM.txt (test_for_verifyBAM.txt is the tester) 


