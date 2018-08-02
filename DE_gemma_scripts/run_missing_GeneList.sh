#!/bin/bash -x

# Before running, make sure you've run make_dirs_byPC_run.pbs to create the directory structure for the by chromosome jobs
# To run, enter './run_short_GeneLists.sh'
	
	GENE_LIST="TreatbySexMainEffect_missing.txt"
	GENE_LIST_DIR=/group/ober-resources/users/mstein3/rna.seq.2017/summary_DE_results/MissingGeneLists
	echo $GENE_LIST_DIR/$GENE_LIST
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_scripts/Run_bimbam_gemma_DE.pbs
	# qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_scripts/Run_bimbam_gemma_MainEffectDE.pbs
	qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_scripts/Run_bimbam_gemma_IntAsMainEffect.pbs