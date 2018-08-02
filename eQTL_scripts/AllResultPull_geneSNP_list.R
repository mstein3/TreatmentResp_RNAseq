#!/usr/bin/env Rscript

# Purpose: From a gene-SNP pair file, pull the gemma output line from that SNP-gene pair. 
#
# Usage: AllResultPull_geneSNP_list.R <gene> <gene_SNP_pair_file> <dir_gemma_covar> <gemma_prefix>
# Submitted with the RunAllResultPull.pbs, using the run_short_Gene_list.sh script. 
# 
# gene: Name of the gene of interest
# 
# gene_SNP_pair_file: File with gene names and desired SNP
#
# dir_gemma_covar: file path to the specific output dir in the gemma directory



### Arguments -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)

f_gene <- args[1]
f_geneSNP_pair <- args[2] 
dir_gemma_covar <- args[3]

#f_gene <- "A2M"
#f_geneSNP_pair <- "/group/ober-resources/users/mstein3/rna.seq.2017/summary_eQTL_results/TMMBatchCor_AgeOnlyCovar/GeneSNP_List_600K_Mash_TMMBatch.txt"
#dir_gemma_covar <- "/scratch/mstein3/eQTL_Hutterites/Male_Null/GemmaDir/TMMbatch"

### Read in the files ------------------------------------

geneSNP_pair <- read.table(f_geneSNP_pair, stringsAsFactors = FALSE, header=FALSE)

f_gemma_output <- paste0("covar1-", f_gene, ".assoc.txt", sep="")

setwd(dir_gemma_covar)
gemma_output <- read.delim(f_gemma_output, stringsAsFactors = FALSE, header=TRUE) 


### Grab the correct line in the data: ------------------------------------

# Find the SNPs to pull out: 
gene_target <- geneSNP_pair[geneSNP_pair$V2 == f_gene,]

rs_to_get<-gene_target[,3]

row_match <- match(rs_to_get, gemma_output$rs)
lines_to_keep <- gemma_output[row_match,]

allRes <- data.frame(f_gene, lines_to_keep, stringsAsFactors = FALSE)

outfile <- sub(".assoc.txt", ".allRes.txt", f_gemma_output)

setwd(dir_gemma_covar)
write.table(allRes, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)











