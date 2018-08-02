#!/usr/bin/env Rscript

# Purpose: Create the gene_SNP list that I'll use to subset from the 32 million to 3.2 million Gene_snp pairs. 

# Usage: CreateSmallerGeneSNPList.R


## TMM Batch ----------------------
# Read in the files needed-------
setwd("/group/ober-resources/users/mstein3/rna.seq.2017/summary_eQTL_results/TMMBatchCor_AgeOnlyCovar")

message("Read in GeneSNP list file")
all_genSNP<-read.table("MN_TMMBatch_GeneSNPList.txt", sep="\t", header=TRUE)

message("Read in top Z score file") 
topres <- read.table("GeneSNPList_LargestZ_TMMBatch.txt", sep="\t", header=FALSE)
colnames(topres) <- c("Gene", "SNP")
topres$GeneSNP_Pair <- paste(topres$Gene, topres$SNP, sep="_")

message("Reorder columns in topres")
# Reorder columns in topres: 
topres2 <- cbind.data.frame(topres$GeneSNP_Pair, topres$Gene, topres$SNP)
colnames(topres2) <- c("GeneSNP_Pair", "Gene", "SNP")

# Organize new gene snp file --------
# Select every 50th SNP: 
small_geneSNP <- all_genSNP[seq(1, nrow(all_genSNP), 50), ]

# Combine with the topres snplist: 
small_geneSNP2 <- rbind.data.frame(small_geneSNP, topres2)

# Remove duplicates: 
small_geneSNP3 <- small_geneSNP2[!duplicated(small_geneSNP2), ]

message("Write out SNPlist file") 
write.table(small_geneSNP3, "GeneSNP_List_600K_Mash_TMMBatch.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
