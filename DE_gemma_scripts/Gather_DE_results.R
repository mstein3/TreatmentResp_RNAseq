#!/usr/bin/env Rscript

# Purpose:Gather the results from each DE .assoc.txt file, and combine into a single file . 

# Usage: Rscript Gather_DE_results.R <dir_gemma_covar> <dir_gemma> <f_top>
# f_top: prefix of output name 

# Arguments -------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)

# dir_gemma_covar
dir_gemma_covar <- args[1]

# dir_gemma
dir_gemma <- args[2]

# top P val file name prefix  
f_top <- args[3]


# Set up combined output file --------------
f_top_full <- file.path(dir_gemma, paste0("CombinedResults_", f_top, ".txt"))


f_top_colnames <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0",
                      "af", "beta", "se", "log1_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score")
cat(f_top_colnames, file = f_top_full, sep = "\t")
cat("\n", file = f_top_full, sep = "", append = TRUE)   


# Functions ----------------
parse_gemma <- function(f_assoc, f_gene, outfile = "") {
  gemma <- read.delim(f_assoc, stringsAsFactors = FALSE)
  top <- data.frame(f_gene, gemma[which.min(abs(gemma$p_lrt)), ],
                    stringsAsFactors = FALSE)
  write.table(top, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  return(invisible(outfile)) 
} 


# Parsing results ----------------------------
setwd(dir_gemma_covar)
file_list <- list.files(pattern="\\.assoc.txt$")
message("Pull result")

for (file in file_list) {
        # Pull out the gene name: 
        # split_by_dash <- strsplit(file, split="NLint_bb-")
        # split_by_dash <- strsplit(file, split="Treat_bb-")
        split_by_dash <- strsplit(file, split="TreatbySex_bb-")
        mat_dash  <- matrix(unlist(split_by_dash), ncol=2, byrow=TRUE)
        gene_assoc <- mat_dash[,2]
        split_by_per <- unlist(strsplit(gene_assoc, split=".assoc"))
        remove_words<-c(".txt")
        gene<- split_by_per[!split_by_per %in% remove_words]
        # Parse file to find top beta
        f_gemma_pval_top <- sub(".assoc.txt", ".withGeneName.txt", file)
        parse_gemma(f_assoc= file, f_gene=gene, outfile = f_gemma_pval_top) 
        }



# Create a combined topPval.txt file: 
message("Combine results")

file_list_2 <- list.files(pattern="\\.withGeneName.txt$")
for (file in file_list_2) {
        f_gemma_pval_top <- file
        # Write to global results file
        system(sprintf("cat %s | sed -e '1d' >> %s", f_gemma_pval_top, f_top_full))
        }





