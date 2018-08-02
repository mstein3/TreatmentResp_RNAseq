#!/usr/bin/env Rscript

# Purpose: Grab and write out the top result line from the topRes files in the gemma dir. Use after MaxResultPull.R, which is generated from the parse_gemma_zscore.R
#
# Usage: Rscript parse_gemma_results.R <dir_gemma> <dir_gemma_covar>

# Arguments -------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)

# dir_gemma_covar
dir_gemma <- args[1]

# top result file (full path)  
dir_gemma_covar <- args[2]

# method used to generate results
method <- args[3]


### Create a global output file --------------

f_top <- file.path(dir_gemma, paste0("topRes-", method, ".txt"))
f_top_colnames <- c("f_gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score")
cat(f_top_colnames, file = f_top, sep = "\t")
cat("\n", file = f_top, sep = "", append = TRUE)    


### Fill the global output file --------------

setwd(dir_gemma_covar)
file_list <- list.files(pattern="\\.topRes.txt$")
message("Pull topRes result")

for (file in file_list) {
        f_gemma_topZ <- file
        # Write to global results file
        system(sprintf("cat %s | sed -e '1d' >> %s", f_gemma_topZ, f_top))
        }

message("Finished the loop")