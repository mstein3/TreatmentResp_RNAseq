# 0) Install and load libraries ====

#call libraries
#library(plyr)
#library(dplyr)
#library(tidyr)
#library(stringr)
#library(ggplot2)
#library(magrittr)
#library(scales)
#library(DescTools)


# 1) Read in files and match them by ids ====

#Set input directory
inp_dir <- "/scratch/mconery/sex-specific/processed_data/"
filter_dir <- "/scratch/mconery/sex-specific/filters/"

#Read in pheno files
pheno_raw_f <- read.table(paste(inp_dir, "females", ".Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pheno_raw_m <- read.table(paste(inp_dir, "males", ".Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pheno_raw <- rbind.data.frame(pheno_raw_f[,1:8], pheno_raw_m[,1:8])
#Read in filter files and filter out FINDIVs not in either ids list
ids_filter <- read.table(paste(filter_dir, "overall.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
pheno_raw <- pheno_raw[which(pheno_raw$Extract_Ord %in% ids_filter$V2),]
remove(ids_filter)
#Sort pheno_file
pheno_raw <- pheno_raw[order(pheno_raw$Extract_Ord),]


#Read in expr file
expr_raw <- read.table(paste(inp_dir, "BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE)
#filter expr files
#Use sprintf to make all numbers four digits, i.e. add leading zeroes if necessary
expr_raw <- expr_raw[,paste("X", pheno_raw$Extract_Ord, sep="")]

#write expression files
for (i in 1:nrow(expr_raw)) {
  #create dfs
  pheno <- as.data.frame(cbind(rep("Hutterites", nrow(pheno_raw)), pheno_raw[,"Extract_Ord"], t(expr_raw[i,])))
  colnames(pheno) <- c("FID", "IID", rownames(expr_raw)[i])
  
  #write files
  write.table(pheno, file = paste("/scratch/mconery/sex-specific/expression_data/", "overall.", rownames(expr_raw)[i], ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
} 

#Write interaction file
write.table(pheno_raw$sex, file = paste(inp_dir, "overall.sex.info.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)
