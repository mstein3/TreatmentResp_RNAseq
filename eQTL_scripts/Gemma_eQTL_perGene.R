#!/usr/bin/env Rscript

# Purpose: Run GEMMA for eQTL analysis on an per gene basis. 

# Usage: Rscript run_gemma.R <counts> <tss> <covar_file> <relat> <prefix_plink> <window> <gene_name> <dir_pheno> <dir_plink> <dir_gemma> <covar_subdir>

# counts: gene-by-sample matrix of normalzied gene expression values
#
# tss: tab-separated file with transcription start site info for each gene in
# counts. No column names, but has the following "chr", "start_tss", "end_gene", "ensg", "coding_type", "gene"
#
# covar: tab-separated file with  covariates.
# First column is 1 for intercept. No column names.
#
# relat: pairwise relatedness values in GEMMA format "id1 id2 value"
#
# prefix_plink: The prefix of the PLINK .bed, .bim, and .fam files.
#
# window: The window to search upstream and downstream of the TSS for SNPs.
#
# gene_name: Name of the gene 
#
# dir_pheno: Top-level directory to save phenotype file for each gene (required
# to subset PLINK data).
#
# dir_plink: Top-level directory to save PLINK files for each gene.
#
# dir_gemma: Top-level directory to store GEMMA results for each gene.
# 
# covar_subdir: name of subdirectory of analysis, useful when running different data sets. Either: covar1, batch_covar1, or sv_covar1


# Arguments----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 11)

# counts
f_exp <- args[1]
# tss
f_tss <- args[2]
# covar
f_covar <- args[3]
# relat
f_relat <- args[4]
# prefix_plink
plink_prefix_all <- args[5]
# window
window <- as.numeric(args[6])
# gene_name
gene_name<-args[7]
# dir_pheno
dir_pheno <- args[8]
# dir_plink
dir_plink <- args[9]
# dir_gemma
dir_gemma <- args[10]
# covar_subdir
covar_subdir <- args[11]

# Add file types to plink prefix:
f_plink_bed <- paste0(plink_prefix_all, ".bed")
f_plink_bim <- paste0(plink_prefix_all, ".bim")
f_plink_fam <- paste0(plink_prefix_all, ".fam")

# Check input ------------------------------------------------------------------

stopifnot(file.exists(f_exp, f_tss, f_covar, f_relat,
                      f_plink_bed, f_plink_bim, f_plink_fam))
                      
# Import data files
exp <- read.delim(f_exp, check.names = FALSE)
tss <- read.table(f_tss, stringsAsFactors = FALSE)
covar <- read.table(f_covar)
# relat <- read.delim(f_relat, stringsAsFactors = FALSE, header = FALSE)
plink_bim <- read.table(f_plink_bim, stringsAsFactors = FALSE)
plink_fam <- read.table(f_plink_fam, stringsAsFactors = FALSE)
                     
stopifnot(nrow(exp) == nrow(tss),
          ncol(exp) == nrow(covar),
          rownames(exp) == tss[, 6])
          
colnames(tss) <- c("chr", "start_tss", "end_gene", "ensg", "coding_type", "gene")         
tss$chr<-as.character(tss$chr) 

if (ncol(covar) == 1) {
  colnames(covar) <- "intercept"
} else if (ncol(covar) == 2) {
		colnames(covar) <- c("intercept", "age")
		} else {
  colnames(covar) <- c("intercept", "age", paste0("covar", seq_along(colnames(covar)[-1:-2])))
} 

rownames(covar) <- colnames(exp)
n_covar <- ncol(covar) - 1

# Check that samples are in PLINK fam file and relatedness matrix
# stopifnot(all(colnames(exp) %in% relat$V1), all(colnames(exp) %in% relat$V1), all(colnames(exp) %in% plink_fam[, 2]))
          
# Make output subdirectories in a separate file.
covar_subdir2 <- paste0(covar_subdir, "/")
dir_pheno_covar <- file.path(dir_pheno, covar_subdir2)
dir_plink_covar <- file.path(dir_plink, covar_subdir2)
dir_gemma_covar <- file.path(dir_gemma, covar_subdir2)

# Create global output file (now done in a separate script)
#f_top <- file.path(dir_gemma, paste0("top-sv-covar-", n_covar, ".txt"))
#f_top_colnames <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0",
#                      "af", "beta", "se", "log1_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score", "n_snps", "n_covar")
#cat(f_top_colnames, file = f_top, sep = "\t")
#cat("\n", file = f_top, sep = "", append = TRUE)    

# Functions for eQTL analysis --------------------------------------------------

get_tss_gene <- function(gene_name, window) {
  stopifnot(is.character(gene_name), is.numeric(window))
  gene_row <- which(tss[,6] == gene_name)
  start_tss <- tss[gene_row, 2]
  chr <- tss[gene_row, 1]
  start <- start_tss - window
  if (start < 0)
    start <- 0
  end = start_tss + window
  out <- c(chr, start, end)
  names(out) <- c("chr", "start", "end")
  return(out)
}  

# Pheno file Format is:
# Column 1/2 = FID/IID
# Column 3 = Phenotype (i.e. gene expression levels)
create_pheno_file <- function(fid, iid, e, f = "") {
  stopifnot(length(fid) == 1 | length(fid) == length(e),
            length(iid) == 1 | length(iid) == length(e),
            is.numeric(e))
  d <- data.frame(fid, iid, e)
  write.table(d, file = f, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  return(invisible(f))
}

create_plink_g <- function(prefix_in, f_pheno, chr, from, to, prefix_out) {
  cmd <- sprintf("plink --bfile %s --make-bed --pheno %s --chr %s --from-bp %s --to-bp  %s --out %s",
                 prefix_in, f_pheno, chr, from, to, prefix_out)
  system(cmd)
  return(invisible(cmd))
}

run_gemma <- function(plink, relatedness, covar, out_prefix, outdir) {
  cmd <- sprintf("gemma -bfile %s -k %s -km 2 -maf 0 -miss 0.2 -c %s -lmm 4 -o %s",
                 plink, relatedness, covar, out_prefix)
  system(cmd)
  # Move to output directory
  cmd2 <- sprintf("mv output/%s* %s", out_prefix, outdir)
  system(cmd2)
  outfile <- paste0(outdir, out_prefix, ".assoc.txt")
  return(invisible(outfile))
}

parse_gemma <- function(f, gene, n_covar, outfile = "") {
  gemma <- read.delim(f, stringsAsFactors = FALSE)
  n_snps <- nrow(gemma)
  top <- data.frame(gene, gemma[which.min(gemma$p_lrt), ], n_snps, n_covar,
                    stringsAsFactors = FALSE)
  write.table(top, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  return(invisible(outfile))
}

# Per gene eQTL analysis -------------------------------------------------------


g <- gene_name
gene_row_num<-which(tss[,6] == gene_name)
iid <- colnames(exp)
fid <- "HUTTERITES"

message("\n\n####\ngene: ", g, "\n####\n\n")

tss_g <- get_tss_gene(gene_name = g, window = window)

f_pheno_g <- paste0(dir_pheno_covar, g, ".pheno")

create_pheno_file(fid = fid, iid = iid, e = as.numeric(exp[gene_row_num, ]),
                    f = f_pheno_g)
                    
plink_prefix_g <- paste0(dir_plink_covar, g)
message("Running PLINK")
 
create_plink_g(prefix_in = plink_prefix_all, f_pheno = f_pheno_g, chr = tss_g["chr"], from = tss_g["start"], to = tss_g["end"], prefix_out = plink_prefix_g)  

if (!file.exists(paste0(plink_prefix_g, ".bed"))) 
    message("\nskipped:\tNo PLINK output\n") 
                 
message("Running GEMMA")
  f_gemma <- run_gemma(plink = plink_prefix_g, relatedness = f_relat, covar = f_covar,
                       out_prefix = paste0("covar", n_covar, "-", g),
                       outdir = dir_gemma_covar)
                       
message("Parse GEMMA")

f_gemma_top <- sub(".assoc.txt", ".top.txt", f_gemma)
parse_gemma(f = f_gemma, gene = g, n_covar = n_covar, outfile = f_gemma_top)

# Write to global results file
# system(sprintf("cat %s | sed -e '1d' >> %s", f_gemma_top, f_top))
          