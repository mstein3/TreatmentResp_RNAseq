# 0) Install and load libraries ====

#install biomaRt
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

#call libraries
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(magrittr)
library(scales)
library(biomaRt)


# 1) Read in pheno file and write ####

#Set input directory
inp_dir <- "/scratch/mconery/sex-specific/processed_data/"

#Read in expr file
x_expr_raw <- read.table(paste(inp_dir, "females.logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = ""), sep = "\t", header = TRUE)
auto_expr_raw <-  read.table(paste(inp_dir, "BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE)

#Get Gene list
genes <- append(rownames(x_expr_raw), rownames(auto_expr_raw))
ensembl=useMart("ensembl")
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","chromosome_name", "start_position"),values=genes,mart= mart)

#Scrub some alternate location entries from G_list
export_list <- G_list[which(substr(G_list$chromosome_name, 1, 3) != "CHR"),]
export_list <- export_list[which(duplicated(export_list$hgnc_symbol == FALSE)),]

#Write list to file
write.table(export_list, file = "/scratch/mconery/sex-specific/filters/gene_locs.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Make genes list files for eQTL loops
x_genes <- export_list[which(export_list$chromosome_name == "X"),"hgnc_symbol"]
auto_genes <- export_list[which(export_list$chromosome_name != "X"),"hgnc_symbol"]
#Insert trailing line feeds (needed to be read in Unix)
x_genes <- append(x_genes, " ")
auto_genes <- append(auto_genes, " ")

#Write genes lists to files
write.table(x_genes, file = "/scratch/mconery/sex-specific/filters/x_genes_eQTL.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(auto_genes, file = "/scratch/mconery/sex-specific/filters/auto_genes_eQTL.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)