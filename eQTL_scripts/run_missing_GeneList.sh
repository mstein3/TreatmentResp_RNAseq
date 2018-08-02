#!/bin/bash -x

# Before running, make sure you've run make_dirs_byPC_run.pbs to create the directory structure for the by chromosome jobs
# To run, enter './run_short_GeneLists.sh'


	
	#GENE_LIST="Uniq_List_UnRunGenes.txt"
	GENE_LIST="TMMBatch_FN_missing.txt"
	GENE_LIST_DIR=/group/ober-resources/users/mstein3/rna.seq.2017/summary_eQTL_results/MissingGeneList
	echo $GENE_LIST_DIR/$GENE_LIST
	 qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_MN.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_ML.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/Run_Gemma_ShortGeneList_All.pbs
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunMaxResultPull.pbs 


	 qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_FN.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_FL.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/Run_Gemma_ShortGeneList_All.pbs

