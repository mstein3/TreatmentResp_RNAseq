#!/usr/bin/env Rscript

# Purpose: Just like pulling the smallest p-value, pull the largest absolute value beta value from each data set: 

# Usage: parse_gemma_betas.R <dir_gemma_covar> <dir_gemma> <analysis method>
# Run with .pbs script: Run_assemble_top_results.pbs

# analysis_method: name of analysis method you'd like to include in the output file title

# Arguments -------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)

dir_gemma_covar <- args[1]
dir_gemma <- args[2] 
method <- args[3]

# To test: 
# dir_gemma_covar <- "/scratch/mstein3/eQTL_Hutterites/Male_Null/GemmaDir/covar1"
# dir_gemma <- "/scratch/mstein3/eQTL_Hutterites/Male_Null/GemmaDir"
# method <- "uncorr"

# Set up top Beta gemma file -------------------
f_top <- file.path(dir_gemma, paste0("top-zscore-", method, ".txt"))


f_top_colnames <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0",
                      "af", "beta", "se", "log1_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score","zscore", "n_snps")
cat(f_top_colnames, file = f_top, sep = "\t")
cat("\n", file = f_top, sep = "", append = TRUE)    



# Functions ----------------
parse_gemma_zscore <- function(f_assoc, f_gene, outfile = "") {
  gemma <- read.delim(f_assoc, stringsAsFactors = FALSE)
  n_snps <- nrow(gemma)
  gemma$zscore <- abs((gemma$beta)/(gemma$se)) 
  top <- data.frame(f_gene, gemma[which.max(abs(gemma$zscore)), ], n_snps,
                    stringsAsFactors = FALSE)
  write.table(top, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  return(invisible(outfile)) 
}




parse_gemma_beta <- function(f_assoc, f_gene, outfile = "") {
  gemma <- read.delim(f_assoc, stringsAsFactors = FALSE)
  n_snps <- nrow(gemma)
  top <- data.frame(f_gene, gemma[which.max(abs(gemma$beta)), ], n_snps,
                    stringsAsFactors = FALSE)
  write.table(top, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  return(invisible(outfile)) 
}



# Parsing results --------------
setwd(dir_gemma_covar)
file_list <- list.files(pattern="\\.assoc.txt$")
message("Pull top zscore result")


# Create the topZ.txt file for each gene: 

for (file in file_list) {
	# Pull out the gene name: 
	split_by_dash <- strsplit(file, split="covar1-")
	mat_dash  <- matrix(unlist(split_by_dash), ncol=2, byrow=TRUE)
	gene_assoc <- mat_dash[,2]
	split_by_per <- unlist(strsplit(gene_assoc, split=".assoc"))
	remove_words<-c(".txt")
	gene<- split_by_per[!split_by_per %in% remove_words]
	# Parse file to find top beta
	f_gemma_beta_top <- sub(".assoc.txt", ".topZ.txt", file)
	parse_gemma_zscore(f_assoc= file, f_gene=gene, outfile = f_gemma_beta_top) 
	}
	


# Create a combined topZ.txt file: 

file_list_2 <- list.files(pattern="\\.topZ.txt$")
for (file in file_list_2) {
        f_gemma_beta_top <- file
        # Write to global results file
        system(sprintf("cat %s | sed -e '1d' >> %s", f_gemma_beta_top, f_top))
        }





