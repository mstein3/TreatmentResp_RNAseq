
#!/usr/bin/env Rscript

# Purpose: Run GEMMA for differential expression by treatment, with a gene by environment (sex) interaction term, using BIMBAM file format. 

# Usage: Usage: Gemma_bimbam_DE.R <counts> <covar_file> <relat> <geno_file> <gxe_file> <gene_name> <dir_pheno> <dir_gemma> 

# counts: gene-by-sample matrix of normalzied gene expression values. Genes as rows, individuals as columns. 

# covar: tab-separated or csv file with  covariates.
# First column is 1 for intercept. No column names.
#
# relat: NbyN matrix of individuals, without row or column names
#
# geno_file: mean genotype file in BIMBAM file format. In this case, the genotype is a treatment. 
#
# gxe_file: column of environmental variable (in this case, sex) 
#
# gene_name: Name of the gene 
#
# dir_pheno: Top-level directory to save phenotype file for each gene 
#
# dir_gemma: Top-level directory to store GEMMA results for each gene.
#
# The main output file is saved as dir_gemma/top-[analysisname]-{name}.txt, where analysis name is
# some attribute of the analysis (number of PCs, test factor used, etc.). 
# gene	chr	rs	ps	n_miss	allele1	allele0	af	beta	se	logl_H1	l_remle	l_mle	p_wald	p_lrt	p_score	n_snps	name

# Arguments----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 8)

# counts
f_exp <- args[1]
# covar
f_covar <- args[2]
# relat
f_relat <- args[3]
# geno_file
f_genotype <- args[4]
# gxe_file
f_gxe <- args[5]
# gene_name
gene_name<-args[6]

# dir_pheno
dir_pheno <- args[7]
# dir_gemma
dir_gemma <- args[8]

#f_exp <- "/group/nicolae-lab/users/mstein3/FcReceptor_June2017/gemma_DE_inputs/input_2geno_Fc/TreatbyGeno89_Model_PoolCorr_GXBIMBAM_51718.txt"
#f_covar <- "/group/nicolae-lab/users/mstein3/FcReceptor_June2017/gemma_DE_inputs/input_2geno_Fc/FcDE_CovarBB_Sex_Age_89.txt"
#f_genotype <- "/group/nicolae-lab/users/mstein3/FcReceptor_June2017/gemma_DE_inputs/input_2geno_Fc/FcDE_rs1801274_bimbam_geno_89.csv"
#gene_name <- "AGRN"
#dir_pheno <- "/group/nicolae-lab/users/mstein3/FcReceptor_June2017/gemma_DE_output/PhenoDir"
#dir_gemma <- "/group/nicolae-lab/users/mstein3/FcReceptor_June2017/gemma_DE_output/GemmaDir"
#f_relat <- "/group/nicolae-lab/users/mstein3/FcReceptor_June2017/additive.coef.2015-12-15"

# Check input ------------------------------------------------------------------

stopifnot(file.exists(f_exp, f_covar, f_relat,
                      f_genotype, f_gxe))
                      
# Import data files

exp0 <- read.delim(f_exp, check.names = FALSE, header=FALSE)
covar <- read.delim(f_covar, header=FALSE)
# genotype <- read.table(f_genotype, stringsAsFactors = FALSE)

gene_names_exp0 <- exp0[,1]

exp <- exp0[,-1]

rownames(exp) <- gene_names_exp0

# exp <- data.frame(exp0[,-1], row.names=exp0[,1])

stopifnot(ncol(exp) == nrow(covar))
# colnames(covar)<-c("intercept", "age", "sex")


# Make output subdirectories in a separate file.
nl_int_subdir <- paste0("TMMBatch_Int", "/")
dir_pheno_covar <- file.path(dir_pheno, nl_int_subdir)
dir_gemma_covar <- file.path(dir_gemma, nl_int_subdir)

# Create global output file
#f_top <- file.path(dir_gemma, paste0("top-bb-nl-int.txt"))
#f_top_colnames <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0",
#                      "af", "beta", "se", "log1_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score", "n_snps")
#cat(f_top_colnames, file = f_top, sep = "\t")
#cat("\n", file = f_top, sep = "", append = TRUE)  

# Functions for differential expression with a genotype analysis --------------------------------------------------

# Pheno file format for bimbam is a column with length = number of samples
# So, just need to pull the row from the gene of interest. 

create_pheno_file <- function(e, f = "") {
  stopifnot(is.numeric(e))
  d <- data.frame(e)
  write.table(d, file = f, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  return(invisible(f))
}


run_gemma <- function(genotype_file, pheno_file, relatedness, covar, gxe, out_prefix, outdir) {
  cmd <- sprintf("gemma -g %s -p %s -k %s -km 1 -notsnp -c %s -gxe %s -lmm 4 -o %s",
                 genotype_file, pheno_file, relatedness, covar, gxe, out_prefix)
  system(cmd)
  # Move to output directory
  cmd2 <- sprintf("mv output/%s* %s", out_prefix, outdir)
  system(cmd2)
  outfile <- paste0(outdir, out_prefix, ".assoc.txt")
  return(invisible(outfile))
}

# Running a bimbam file for gemma: 
# gemma -g [mean genotype file] -p [phenotype file name] -k [relatedness matrix] -c [covar] -gxe [gxe_file] -lmm 4 -notsnp -o [output] 


# Per gene DE analysis -------------------------------------------------------

g <- gene_name
gene_row_num<-which(rownames(exp) == gene_name)

message("\n\n####\ngene: ", g, "\n####\n\n")

f_pheno_g <- paste0(dir_pheno_covar, g, ".pheno")

create_pheno_file(e = as.numeric(exp[gene_row_num,]),
                    f = f_pheno_g)
                    
message("Running GEMMA")
  f_gemma <- run_gemma(genotype_file = f_genotype, pheno_file = f_pheno_g, relatedness = f_relat, covar = f_covar, gxe = f_gxe,
                       out_prefix = paste0("NLint_bb-", g),
                       outdir = dir_gemma_covar)
                       



          




