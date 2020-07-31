# 0) Install and load libraries ====

#install limma 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("limma")

#call libraries
library(limma)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(magrittr)
library(scales)
library(ggrepel)
library(edgeR)


# 1) Read in Required STAR Files ====

#Set input directory
inp_dir <- "C:/Users/Janel/Documents/sex-specific/raw_data/"
process_dir <- "C:/Users/Janel/Documents/sex-specific/processed_data/"
filter_dir <- "C:/Users/Janel/Documents/sex-specific/filters/"

#Read in STAR gene counts for all genes
STAR_gene_counts <- read.table(file = paste(inp_dir, "RawGeneCounts_RunsCombined_ENSG_names_noQC_9117.txt", sep = ""), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
#Read in ENSG key
ENSG_key <- read.table(file = paste(inp_dir, "ENSG_GeneNames_Chr_All.txt", sep = ""), sep = " ", header = FALSE, stringsAsFactors = FALSE)
#Read in list of processed data to get final list of qc'ed genes and kept extract_ords
expr_original <- read.table(file = paste(process_dir, "BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

#Filter ENSG_key for genes that pass final qc filtersn (old version)
#ENSG_key <- ENSG_key[which(ENSG_key$V6 %in% rownames(expr_original)),]

#Convert STAR_gene_counts row names using the ENSG key
#Rewrite rownames of ENSG_key with ENSG ids
rownames(ENSG_key) <- ENSG_key[,"V4"]

#Create cycle list for apply function
cycle_list <- ENSG_key[which(duplicated(ENSG_key$V6) == FALSE),"V6"]

#Filter STAR_gene_counts for extract ords that pass qc filters
STAR_gene_counts <- STAR_gene_counts[,sprintf("X%04d", as.numeric(substr(colnames(expr_original), 2, nchar(colnames(expr_original)))))]

#Create function to sum rows by taking row name of a row in the processed file 
sum_raw_rows <- function(gene_symbol, ENSG_key, STAR_gene_counts){
  
  #Filter ENSG_key
  ENSG_key <- ENSG_key[which(ENSG_key$V6 == gene_symbol),]
  
  #Filter STAR_gene_counts for ENSGs of interest
  STAR_gene_counts <- STAR_gene_counts[rownames(ENSG_key),]
  
  #Temp store variable
  temp_hold <- colSums(STAR_gene_counts)
  
  #Return vector of column sums
  return(temp_hold)
  
}

#apply sum_raw_rows on vector of row names of STAR_gene_counts_processed to create processed_data_frame
STAR_gene_counts_processed <- as.data.frame(t(sapply(cycle_list, sum_raw_rows, ENSG_key=ENSG_key, STAR_gene_counts=STAR_gene_counts)))

#Clean-up junk 
suppressWarnings(remove(STAR_gene_counts, ENSG_key, sum_raw_rows, cycle_list))


# 2) Check HLA-DRB3, HLA-DRB4, and HLA-DRB5 reads ====

#Not everyone should have reads for these three genes as not everyone has a copy of all three

#Filter for the three genes and ATF1 (Which will be a reference point)
DRB_Gene_Counts <- STAR_gene_counts_processed[which(rownames(STAR_gene_counts_processed) %in% c("HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "ATF1")),]

#It seems like there are a few individuals with almost no mapped reads to the DRB gene. Low genetic diversity in the population and the 
#similarity of the DRB genes makes me think these results are reasonable.  The low reads are also paired by FINDIVs (1153/1104 and 1110/1261),
#so the low read is unaffected by treatment as expected.  We are going to continue to use the HLA-DRB5 reads.  

#Clean-up
remove(DRB_Gene_Counts)

# 3) Read in New HLA data from Selene ====

#Read in Selene's Data
HLA_gene_counts <- read.table(file = paste(inp_dir, "HLA_gene_counts.tsv", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Remove non-extract_ord info in first column
HLA_gene_counts[,"subject"] <- substr(HLA_gene_counts[,"subject"], 1, regexpr("_", HLA_gene_counts[,"subject",]) - 1)

#Honestly this line does nothing
HLA_gene_counts$subject <- as.character(as.integer(HLA_gene_counts$subject))

#Create vector of final extract ords
extract_ords <- as.character(as.integer(substr(colnames(STAR_gene_counts_processed), 2, nchar(STAR_gene_counts_processed))))
#Create vector of HLA gene names
HLA_genes <- HLA_gene_counts$locus[which(duplicated(HLA_gene_counts$locus) == FALSE)]

#Define internal function to combine rows for a given HLA gene for a specific individual
internal_combine_counts <- function(HLA_gene, HLA_gene_counts){
  
  #Filter HLA_gene_counts for HLA_gene of interest
  HLA_gene_counts <- HLA_gene_counts[which(HLA_gene_counts$locus == HLA_gene),]
  
  #Return sum of columns
  return(sum(HLA_gene_counts$est_counts))
  
}

#Define function to combine rows per person in extract_ords vector
combine_counts <- function(extract_ord, HLA_genes, HLA_gene_counts){
  
  #Filter HLA_gene_counts for extract_ord of interest
  HLA_gene_counts <- HLA_gene_counts[which(HLA_gene_counts$subject == extract_ord),]
  
  #Apply internal function to get counts for the extract ord and return vector
  return(sapply(HLA_genes, internal_combine_counts, HLA_gene_counts=HLA_gene_counts))
  
}

#Apply outer function to list of extract_ords
HLA_processed_counts <- as.data.frame(sapply(extract_ords, combine_counts, HLA_genes=HLA_genes, HLA_gene_counts=HLA_gene_counts))

#Clean-up
suppressWarnings(remove(combine_counts, internal_combine_counts, miss_extract_ords, HLA_gene_counts, extract_ords))



## Make bar graphs comparing new and old raw HLA counts
#Extract HLA rows from STAR gene counts
HLA_STAR_counts <- STAR_gene_counts_processed[rownames(HLA_processed_counts),]
#Need to make a single matrix with columns: gene name, extract_ord, Updated vs. STAR
#First create blank data frame
barplot_df <- data.frame(matrix(data = NA, nrow = 0, ncol = 4))
colnames(barplot_df) <- c("gene", "extract_ord", "method", "count")
#Loop through rows in order to create the data frame
for (i in 1:nrow(HLA_STAR_counts)) {
  #Create temp_df to hold data from HLA_STAR_counts
  temp_df <- cbind.data.frame(gene=as.character(rep(rownames(HLA_STAR_counts)[i], ncol(HLA_STAR_counts))),
                              extract_ord=as.character(substr(colnames(HLA_STAR_counts), 2, nchar(colnames(HLA_STAR_counts)))),
                              method=as.character(rep("original", ncol(HLA_STAR_counts))),
                              count=t(HLA_STAR_counts[i,]))
  colnames(temp_df)[4] <- "count"
  #Rbind temp_df to barplot_df
  barplot_df <- rbind.data.frame(barplot_df, temp_df)
  #Repeat process for HLA_processed_counts
  temp_df <- cbind.data.frame(gene=as.character(rep(rownames(HLA_processed_counts)[i], ncol(HLA_processed_counts))),
                              extract_ord=sprintf("%04d", as.numeric(colnames(HLA_processed_counts))),
                              method=as.character(rep("updated", ncol(HLA_processed_counts))),
                              count=t(HLA_processed_counts[i,]))
  colnames(temp_df)[4] <- "count"
  #Rbind temp_df to barplot_df
  barplot_df <- rbind.data.frame(barplot_df, temp_df)
}
remove(temp_df)
#Set out_dir
out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/HLA_barplots/"
#Make plots for all six genes by cycling throuh HLA_genes
for (HLA_gene in HLA_genes) {
  #filter for gene of interest
  barplot_df_filtered <- barplot_df[which(barplot_df$gene == HLA_gene),]
  #make boxplot
  tiff(paste(out_dir, HLA_gene, ".comparison.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  print(ggplot(data=barplot_df_filtered, aes(x=extract_ord, y=count, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()))
  dev.off()
}
remove(barplot_df_filtered, barplot_df, HLA_STAR_counts)

#Filter out extract ords missing HLA 
HLA_processed_counts <- HLA_processed_counts[,which(colSums(HLA_processed_counts) != 0)]

#Overwrite column names of HLA_processed_counts and expr_original
colnames(HLA_processed_counts) <- sprintf("X%04d", as.numeric(colnames(HLA_processed_counts)))
colnames(expr_original) <- sprintf("X%04d", as.numeric(substr(colnames(expr_original), 2, nchar(colnames(expr_original)))))

#Sort and filter columns in three matrices
HLA_processed_counts <- HLA_processed_counts[,order(colnames(HLA_processed_counts))]
expr_original <- expr_original[,colnames(HLA_processed_counts)]
STAR_gene_counts_processed <- STAR_gene_counts_processed[,colnames(HLA_processed_counts)]

#Overwrite rows of STAR_gene_counts_processed with new HLA data
STAR_gene_counts_processed[rownames(HLA_processed_counts),] <- HLA_processed_counts


# 4) Run Kevin's script for normalization ====

#Filter for 12,478 genes in original expression matrix (new code added on second pass)
#Running with Kevin's instructions, we get the following discrepancies:
#208 genes among original 12,478 are no longer in the final counts matrix
#2,000+ genes are now included that were not previously

#Convert to counts per million and filter for genes with cpm>1 in 20 or more samples (per methods)
STAR_gene_counts_final <- STAR_gene_counts_processed[rownames(expr_original),]

#TMM normalization
# To calculate the TMM normalization factors, create a DGElist using the edgeR package: 
dge.nl.nonorm <- DGEList(counts=STAR_gene_counts_final)
# Perform TMM normalization using the calcNormFactors function:
dge.nl <- calcNormFactors(dge.nl.nonorm)

#Voom transformation
# Apply voom transformation
v.nl <- voom(dge.nl, design=NULL,plot=TRUE)
v.gx <- v.nl$E

#Create expr_raw
expr_raw <- as.data.frame(v.gx)

#Clean-up
suppressWarnings(remove(dge.nl, dge.nl.nonorm, STAR_gene_counts_processed, STAR_gene_counts_final, v.gx, v.nl, i))



# 5) Remove covariates from expr_raw ====

#Read in pheno_file
pheno_raw <- read.table(file = paste(inp_dir, "Pheno_Rready_403comb_8117.txt", sep = ""), sep = "\t", header = TRUE, row.names = 2, stringsAsFactors = FALSE)
#Reformat extract_ord rownames
rownames(pheno_raw) <- sprintf("X%04d", as.numeric(rownames(pheno_raw)))
pheno_raw <- pheno_raw[order(rownames(pheno_raw)),]
#Filter pheno raw for extract orders of interest
pheno_raw <- pheno_raw[colnames(expr_raw),]

#Recode extract batch and pool code as factors
pheno_raw$Extract_Batch <- as.factor(pheno_raw$Extract_Batch)
pheno_raw$Pool_1.0 <- as.factor(pheno_raw$Pool_1.0)

#Remove batch effects
expr_processed <- as.data.frame(removeBatchEffect(expr_raw, batch = pheno_raw$Pool_1.0, batch2 = pheno_raw$Extract_Batch))

#clean-up
suppressWarnings(remove(expr_raw, HLA_processed_counts))

# 6) Compare expression matricies and isolate HLA genes for export ====

#Filter expr_processed for genes of interest
expr_processed_overlap <- expr_processed[rownames(expr_original),]

#Sort rows the same way
expr_original <- expr_original[order(rownames(expr_original)),]
expr_processed_overlap <- expr_processed_overlap[order(rownames(expr_processed_overlap)),]

# Calculate percent of values deviating by more than 10%
temp <- rowSums(abs((expr_processed_overlap - expr_original) / expr_original) > 0.10)
comparison_ten <- sum(rowSums(abs((expr_processed_overlap[which(is.na(temp) == FALSE),] - expr_original[which(is.na(temp) == FALSE),]) / expr_original[which(is.na(temp) == FALSE),]) > 0.10)) / (nrow(expr_original[which(is.na(temp) == FALSE),]) * ncol(expr_original[which(is.na(temp) == FALSE),]))
comparison_twenty <- sum(rowSums(abs((expr_processed_overlap[which(is.na(temp) == FALSE),] - expr_original[which(is.na(temp) == FALSE),]) 
                                     / expr_original[which(is.na(temp) == FALSE),]) > 0.20)) / (nrow(expr_original[which(is.na(temp) == FALSE),]) * ncol(expr_original[which(is.na(temp) == FALSE),]))
#Only 14% of expression values vary by more than 10% and only 9% vary by more than 20%.  

#clean-up
suppressWarnings(remove(expr_processed_overlap, comparison_ten, comparison_twenty, temp, comparison))

#Filter expr_processed for 8 HLA genes 
HLA_expr <- expr_processed[HLA_genes,]

#Rename column names
colnames(HLA_expr) <- as.character(as.numeric(substr(colnames(HLA_expr),2,nchar(colnames(HLA_expr)))))

#Write to file
write.table(HLA_expr, file = paste(process_dir, "HLA_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#Create male and female iid lists to write to filter dir
male_iids <- cbind.data.frame(fid=rep("HUTTERITES", nrow(pheno_raw[which(pheno_raw$sex == 0),])), iid=pheno_raw[which(pheno_raw$sex == 0),"FINDIV"])
female_iids <- cbind.data.frame(fid=rep("HUTTERITES", nrow(pheno_raw[which(pheno_raw$sex == 1),])), iid=pheno_raw[which(pheno_raw$sex == 1), "FINDIV"])

#Sort id lists  
male_iids <- male_iids[order(male_iids$iid),]
female_iids <- female_iids[order(female_iids$iid),]

#Remove duplicates from lists
male_iids <- male_iids[which(duplicated(male_iids) == FALSE),]
female_iids <- female_iids[which(duplicated(female_iids) == FALSE),]

#Write to file
write.table(male_iids, file = paste(filter_dir, "HLA.males.iids.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(female_iids, file = paste(filter_dir, "HLA.females.iids.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Write list of HLA_genes to file too
write.table(HLA_genes, file = paste(filter_dir, "HLA_auto_genes_DE.tsv", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
