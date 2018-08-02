#!/bin/bash -x

# to submit this, enter './adaptor_list_through.sh' 
# the read command is going to look at whatever file you pipe into the script with the < operator
while read LINE; do
	echo "$LINE"
	#qsub -v FILE=$LINE fastqc_run.pbs
	qsub -v FILE=$LINE star_run.pbs
	#qsub -v FILE=$LINE multiQC_run.pbs
	#qsub -v FILE=$LINE fastqc_leftovers.pbs
	
done < $1	
#done < /group/ober-resources/users/mstein3/RNAseq/Processing_scripts/adaptor_index_list


# Additional files run outside of this job submission framework (in order): 
# filelength.pbs (run whenever). 
# fastqc_summary.py to summarize fastQC output. (run after Run_fastQC.pbs is completed)
# rmlc.sh to remove gene count files that have no counts/data (after Gene_counts.pbs is completed, required to merge gene counts in R later). 
