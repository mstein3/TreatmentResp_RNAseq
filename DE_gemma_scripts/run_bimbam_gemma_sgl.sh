#!/bin/bash -x

# Before running, make sure you've set up the output directory structure.
# To run, enter './run_short_GeneLists.sh'

 for i in $(seq 1 1 125); 
# for i in $(seq 1 1 2); 
do 
	
	GENE_LIST="Short_gene_list_${i}.txt"
	GENE_LIST_DIR=/group/ober-resources/users/mstein3/rna.seq.2017/input_files_eQTL/Gene_Lists
	echo $GENE_LIST
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_scripts/Run_bimbam_gemma_DE.pbs
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_scripts/Run_bimbam_gemma_MainEffectDE.pbs
	qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_scripts/Run_bimbam_gemma_IntAsMainEffect.pbs
	
done

