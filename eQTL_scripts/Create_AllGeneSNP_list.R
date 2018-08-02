#!/usr/bin/env Rscript

# Purpose:  Gathering all gene-snp pairs tested for each gene from the GEMMA output. Then, put all of those files together into a single file. 
# Based off of Create_AllBetas_AllSE.R, where the geneSNP list was pulled, and the betas, SEs, and p-values. 
# Also based off of Combine_BetaSE_files.R, where betaSE files are combined into one file.  

# Usage: Rscript Create_AllGeneSNP_list.R <dir_gemma_covar> <dir_gemma> <collected_output_name>
# Run in all 4 datasets, in case it's not the same. 

# Arguments -------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)

dir_gemma_covar <- args[1]

dir_gemma <- args[2]

collected_output_name <- args[3]

# Pull beta and se for each gene, with a gene-snp column ------------
setwd(dir_gemma_covar)
file_list <- list.files(pattern="\\.assoc.txt$")

message("GeneSNP data from each assoc file")

for (file in file_list) {
        split_by_dash<-strsplit(file, split="covar1-")
        mat_dash  <- matrix(unlist(split_by_dash), ncol=2, byrow=TRUE)
        gene_assoc <-mat_dash[,2]
        split_by_per <- unlist(strsplit(gene_assoc, split=".assoc"))
        remove_words<-c(".txt")
        gene_name<- split_by_per[!split_by_per %in% remove_words]
        
        assoc <- read.delim(file, stringsAsFactors = FALSE, header=TRUE)
        snp_only <- assoc$rs
        length_gene <- nrow(assoc)
        gene_list <- rep(gene_name, length_gene)
        gene_snp <- paste(gene_list, snp_only, sep="_")
        table1 <- cbind.data.frame(gene_snp, gene_list, snp_only)
                table_name <- paste(gene_name, ".GeneSNP.txt", sep="")
                write.table(table1, table_name, sep="\t", quote=FALSE, row.names=FALSE)
		}


# Create the combined beta se output file --------------

f_betaSE <- file.path(dir_gemma, paste0(collected_output_name,".txt", sep=""))

f_betaSE_colnames <- c("GeneSNP_Pair", "Gene", "SNP")
cat(f_betaSE_colnames, file = f_betaSE, sep = "\t")
cat("\n", file = f_betaSE, sep = "", append = TRUE)    
        
                
# Combine all of the betaSE files -------------         
setwd(dir_gemma_covar)
file_list <- list.files(pattern="\\.GeneSNP.txt$")
message("Combine Gene SNP files")                
                
for (file in file_list) {
        f_file <- file
        # Write to global results file
        system(sprintf("cat %s | sed -e '1d' >> %s", f_file, f_betaSE))
        }                       

message("finished")



