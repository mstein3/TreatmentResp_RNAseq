#!/usr/bin/env Rscript

# Purpose: From a gene-SNP pair file, pull the gemma output line from that SNP-gene pair. 
#
# Usage: MaxResultPull_geneSNP_list.R <gene> <gene_SNP_pair_file> <dir_gemma_covar> <gemma_prefix>
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


### Read in the files ------------------------------------

geneSNP_pair <- read.table(f_geneSNP_pair, stringsAsFactors = FALSE, header=FALSE)
rownames(geneSNP_pair) <- geneSNP_pair[,1]

f_gemma_output <- paste0("covar1-", f_gene, ".assoc.txt", sep="")

setwd(dir_gemma_covar)
gemma_output <- read.delim(f_gemma_output, stringsAsFactors = FALSE, header=TRUE) 


### Grab the correct line in the data: ------------------------------------

# Find the SNP to pull out: 
rs_to_get<-geneSNP_pair[f_gene,2]

row_match <- match(rs_to_get, gemma_output$rs)
line_to_keep <- gemma_output[row_match,]

top_beta <- data.frame(f_gene, line_to_keep, stringsAsFactors = FALSE)

outfile <- sub(".assoc.txt", ".topRes.txt", f_gemma_output)

setwd(dir_gemma_covar)
write.table(top_beta, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)











