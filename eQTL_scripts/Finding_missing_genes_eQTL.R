### Collecting the genes that didn't run in a cis-eQTL analysis with gemma
## Then, make a gene list of genes that didn't run 

# Usage: Rscript Finding_missing_genes_eQTL.R <dir_gemma_covar> <dir_gene_output> <output_prefix>


# Arguments -------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
# dir_gemma_covar
dir_gemma_covar <- args[1]

dir_gene_output <- args[2]

output_prefix <- args[3]


# Testing example
#dir_gemma_covar <- "/scratch/mstein3/eQTL_Hutterites/Female_Null/GemmaDir/sv_covar1"
#dir_gene_output <- "/scratch/mstein3/"
#output_prefix <- "test_delete_later"


# Parsing results ----------------------------
setwd(dir_gemma_covar)
file_list <- list.files(pattern="\\.assoc.txt$")
message("Split_files_into_genes")


split_by_dash<-strsplit(file_list, split="covar1-")
mat_dash  <- matrix(unlist(split_by_dash), ncol=2, byrow=TRUE)
gene_assoc <-mat_dash[,2]
split_by_per <- unlist(strsplit(gene_assoc, split=".assoc"))

remove_words<-c(".txt")

gene_list<- split_by_per[!split_by_per %in% remove_words]
gene_list <- as.data.frame(gene_list)

# Read in full gene list: 
genes_all <- read.table("/group/ober-resources/users/mstein3/rna.seq.2017/input_files_eQTL/GeneList_JobConf.txt")

message("Determine which genes are missing")

# Determine which genes are missing: 
missing_genes <- setdiff(genes_all[,1], gene_list[,1])
missing_genes <- as.data.frame(missing_genes)

# Write out this as a gene file----------------------
setwd(dir_gene_output) 

write.table(missing_genes, paste0(output_prefix, "_missing.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)