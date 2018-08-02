#!/bin/bash -x

# Before running, make sure you've run make_dirs_byPC_run.pbs to create the directory structure for the by chromosome jobs
# To run, enter './run_short_GeneLists.sh'

for i in $(seq 1 1 125); 
#for i in $(seq 1 1 2); 
do 
	
	GENE_LIST="Short_gene_list_${i}.txt"
	GENE_LIST_DIR=/group/ober-resources/users/mstein3/rna.seq.2017/input_files_eQTL/Gene_Lists
	echo $GENE_LIST
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_MN.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_ML.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/Run_Gemma_ShortGeneList_All.pbs
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunMaxResultPull.pbs 
	 qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunAllResultPull.pbs 
	
done

#for i in $(seq 1 1 125); 
#do 
	
	# GENE_LIST="Short_gene_list_${i}.txt"
	# GENE_LIST_DIR=/group/ober-resources/users/mstein3/rna.seq.2017/input_files_eQTL/Gene_Lists
	# echo $GENE_LIST
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_FN.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/RunGemma_ShGeneL_FL.pbs 
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/eQTL_scripts/Run_Gemma_ShortGeneList_All.pbs
	
#done