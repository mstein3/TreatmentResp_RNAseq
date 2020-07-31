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
library(DescTools)
library(ggrepel)
library(ggpubr)
library(VennDiagram)



# 1) Read in results and make Venn Diagrams ====

#Name junk HLA_genes for future use
junk_genes <- c("HLA-DPA1", "HLA-DRB5", "HLA-G")

#Read in file
mashr_raw <- read.table(file = "C:/Users/Janel/Documents/sex-specific/eQTL_analysis/AllRes_lfsr_combined_031120.txt", 
                        stringsAsFactors = FALSE, header = TRUE)

#Add gene column to mashr_raw
mashr_raw <- cbind.data.frame(gene=as.character(rownames(mashr_raw)), mashr_raw)
mashr_raw$gene <- as.character(mashr_raw$gene)
mashr_raw$gene <- substr(mashr_raw$gene,1,regexpr("\\:", mashr_raw$gene) - 1)

#Filter out junk genes
mashr_raw <- mashr_raw[which(mashr_raw$gene %in% junk_genes == FALSE),]

#Read in gene lists
auto_genes <- read.table(file = "C:/Users/Janel/Documents/Sex-specific/processed_data/BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = "\t", stringsAsFactors = FALSE)
auto_genes <- rownames(auto_genes)
x_genes <- read.table(file = "C:/Users/Janel/Documents/Sex-specific/processed_data/females.logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = "\t", stringsAsFactors = FALSE)
x_genes <- rownames(x_genes)

#Filter auto and x results 
auto_eQTLs <- mashr_raw[which(mashr_raw$gene %in% auto_genes),c(2:5)]
x_eQTLs <- mashr_raw[which(mashr_raw$gene %in% x_genes),c(2:5)]

#Count number of significant gene-snp combos in at least one condition
num_eQTLs <- nrow(mashr_raw[which(mashr_raw$MaleNull < 0.05 | mashr_raw$MaleLPS < 0.05 | mashr_raw$FemaleNull < 0.05 | mashr_raw$FemaleLPS < 0.05),])
num_auto_eQTLs <- nrow(auto_eQTLs[which(auto_eQTLs$MaleNull < 0.05 | auto_eQTLs$MaleLPS < 0.05 | auto_eQTLs$FemaleNull < 0.05 | auto_eQTLs$FemaleLPS < 0.05),])
num_x_eQTLs <- nrow(x_eQTLs[which(x_eQTLs$MaleNull < 0.05 | x_eQTLs$MaleLPS < 0.05 | x_eQTLs$FemaleNull < 0.05 | x_eQTLs$FemaleLPS < 0.05),])

#Plot X-linked eQTLs
#Calculate overlap groups
overlap_values<-get.venn.partitions(list("Females/LPS" = rownames(x_eQTLs[which(x_eQTLs$FemaleLPS < 0.05),]), 
                                         "Females/VEH" = rownames(x_eQTLs[which(x_eQTLs$FemaleNull < 0.05),]),
                                         "Males/LPS" = rownames(x_eQTLs[which(x_eQTLs$MaleLPS < 0.05),]), 
                                         "Males/VEH" = rownames(x_eQTLs[which(x_eQTLs$MaleNull < 0.05),])))
#Redefine output diagram for venn diagrams
out_dir <- "C:/Users/Janel/Documents/Sex-specific/Visuals/venn_diagrams/"
#Draw Venn Diagram
tiff(paste(out_dir, "mashr_x_venn_diagram.jpg", sep = ""), width = 11.5, height = 8, units = 'in', res = 500)
print(draw.quad.venn(area1 = length(x_eQTLs$FemaleLPS[which(x_eQTLs$FemaleLPS < 0.05)]), 
                     area2 = length(x_eQTLs$FemaleNull[which(x_eQTLs$FemaleNull < 0.05)]), 
                     area3 = length(x_eQTLs$MaleLPS[which(x_eQTLs$MaleLPS < 0.05)]), 
                     area4 = length(x_eQTLs$MaleNull[which(x_eQTLs$MaleNull < 0.05)]), 
               n12 = overlap_values$..count..[1] + overlap_values$..count..[5] + overlap_values$..count..[9] + overlap_values$..count..[13], 
               n13 = overlap_values$..count..[1] + overlap_values$..count..[3] + overlap_values$..count..[9] + overlap_values$..count..[11], 
               n14 = overlap_values$..count..[1] +  overlap_values$..count..[3] + overlap_values$..count..[5] + overlap_values$..count..[7], 
               n23 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[9] + overlap_values$..count..[10], 
               n24 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[5] + overlap_values$..count..[6], 
               n34 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[3] + overlap_values$..count..[4], 
               n123 = overlap_values$..count..[1] + overlap_values$..count..[9], 
               n124 = overlap_values$..count..[1] + overlap_values$..count..[5],
               n134 = overlap_values$..count..[1] + overlap_values$..count..[3], 
               n234 = overlap_values$..count..[1] + overlap_values$..count..[2], 
               n1234 = overlap_values$..count..[1],
               category = c("Females/LPS", "Females/VEH", "Males/LPS", "Males/VEH"), 
               fill = c("red", "red", "skyblue1", "skyblue1"),
               alpha = c(0.9, 0.2, 0.9, 0.4), cex = rep(2.5, 15), cat.cex = rep(2.5, 4),
               cat.just = list(c(0.4,0), c(0.6,0), c(0.4,0), c(0.6,0))))
dev.off()
#Make PDF for Michelle
pdf(paste(out_dir, "Figure_5B.pdf", sep = ""), width = 11.5, height = 8, useDingbats = FALSE)
print(draw.quad.venn(area1 = length(x_eQTLs$FemaleLPS[which(x_eQTLs$FemaleLPS < 0.05)]), 
                     area2 = length(x_eQTLs$FemaleNull[which(x_eQTLs$FemaleNull < 0.05)]), 
                     area3 = length(x_eQTLs$MaleLPS[which(x_eQTLs$MaleLPS < 0.05)]), 
                     area4 = length(x_eQTLs$MaleNull[which(x_eQTLs$MaleNull < 0.05)]), 
                     n12 = overlap_values$..count..[1] + overlap_values$..count..[5] + overlap_values$..count..[9] + overlap_values$..count..[13], 
                     n13 = overlap_values$..count..[1] + overlap_values$..count..[3] + overlap_values$..count..[9] + overlap_values$..count..[11], 
                     n14 = overlap_values$..count..[1] +  overlap_values$..count..[3] + overlap_values$..count..[5] + overlap_values$..count..[7], 
                     n23 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[9] + overlap_values$..count..[10], 
                     n24 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[5] + overlap_values$..count..[6], 
                     n34 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[3] + overlap_values$..count..[4], 
                     n123 = overlap_values$..count..[1] + overlap_values$..count..[9], 
                     n124 = overlap_values$..count..[1] + overlap_values$..count..[5],
                     n134 = overlap_values$..count..[1] + overlap_values$..count..[3], 
                     n234 = overlap_values$..count..[1] + overlap_values$..count..[2], 
                     n1234 = overlap_values$..count..[1],
                     category = c("Females/LPS", "Females/VEH", "Males/LPS", "Males/VEH"), 
                     fill = c("red", "red", "skyblue1", "skyblue1"),
                     alpha = c(0.9, 0.2, 0.9, 0.4), cex = rep(2.5, 15), cat.cex = rep(2.5, 4),
                     cat.just = list(c(0.4,0), c(0.6,0), c(0.4,0), c(0.6,0))))
dev.off()


#Plot Autosomal eQTLs
#Calculate overlap groups
overlap_values<-get.venn.partitions(list("Females/LPS" = rownames(auto_eQTLs[which(auto_eQTLs$FemaleLPS < 0.05),]), 
                                         "Females/VEH" = rownames(auto_eQTLs[which(auto_eQTLs$FemaleNull < 0.05),]),
                                         "Males/LPS" = rownames(auto_eQTLs[which(auto_eQTLs$MaleLPS < 0.05),]), 
                                         "Males/VEH" = rownames(auto_eQTLs[which(auto_eQTLs$MaleNull < 0.05),])))
#Redefine output diagram for venn diagrams
out_dir <- "C:/Users/Janel/Documents/Sex-specific/Visuals/venn_diagrams/"
#Draw Venn Diagram
tiff(paste(out_dir, "mashr_auto_venn_diagram.jpg", sep = ""), width = 11.5, height = 8, units = 'in', res = 500)
print(draw.quad.venn(area1 = length(auto_eQTLs$FemaleLPS[which(auto_eQTLs$FemaleLPS < 0.05)]), 
                     area2 = length(auto_eQTLs$FemaleNull[which(auto_eQTLs$FemaleNull < 0.05)]), 
                     area3 = length(auto_eQTLs$MaleLPS[which(auto_eQTLs$MaleLPS < 0.05)]), 
                     area4 = length(auto_eQTLs$MaleNull[which(auto_eQTLs$MaleNull < 0.05)]), 
                     n12 = overlap_values$..count..[1] + overlap_values$..count..[5] + overlap_values$..count..[9] + overlap_values$..count..[13], 
                     n13 = overlap_values$..count..[1] + overlap_values$..count..[3] + overlap_values$..count..[9] + overlap_values$..count..[11], 
                     n14 = overlap_values$..count..[1] +  overlap_values$..count..[3] + overlap_values$..count..[5] + overlap_values$..count..[7], 
                     n23 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[9] + overlap_values$..count..[10], 
                     n24 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[5] + overlap_values$..count..[6], 
                     n34 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[3] + overlap_values$..count..[4], 
                     n123 = overlap_values$..count..[1] + overlap_values$..count..[9], 
                     n124 = overlap_values$..count..[1] + overlap_values$..count..[5],
                     n134 = overlap_values$..count..[1] + overlap_values$..count..[3], 
                     n234 = overlap_values$..count..[1] + overlap_values$..count..[2], 
                     n1234 = overlap_values$..count..[1],
                     category = c("Females/LPS", "Females/VEH", "Males/LPS", "Males/VEH"), 
                     fill = c("red", "red", "skyblue1", "skyblue1"),
                     alpha = c(0.9, 0.2, 0.9, 0.4), cex = rep(2.5, 15), cat.cex = rep(2.5, 4),
                     cat.just = list(c(0.4,0), c(0.6,0), c(0.4,0), c(0.6,0))))
dev.off()
#Make PDF for Michelle
pdf(paste(out_dir, "Figure_5A.pdf", sep = ""), width = 11.5, height = 8, useDingbats = FALSE)
print(draw.quad.venn(area1 = length(auto_eQTLs$FemaleLPS[which(auto_eQTLs$FemaleLPS < 0.05)]), 
                     area2 = length(auto_eQTLs$FemaleNull[which(auto_eQTLs$FemaleNull < 0.05)]), 
                     area3 = length(auto_eQTLs$MaleLPS[which(auto_eQTLs$MaleLPS < 0.05)]), 
                     area4 = length(auto_eQTLs$MaleNull[which(auto_eQTLs$MaleNull < 0.05)]), 
                     n12 = overlap_values$..count..[1] + overlap_values$..count..[5] + overlap_values$..count..[9] + overlap_values$..count..[13], 
                     n13 = overlap_values$..count..[1] + overlap_values$..count..[3] + overlap_values$..count..[9] + overlap_values$..count..[11], 
                     n14 = overlap_values$..count..[1] +  overlap_values$..count..[3] + overlap_values$..count..[5] + overlap_values$..count..[7], 
                     n23 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[9] + overlap_values$..count..[10], 
                     n24 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[5] + overlap_values$..count..[6], 
                     n34 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[3] + overlap_values$..count..[4], 
                     n123 = overlap_values$..count..[1] + overlap_values$..count..[9], 
                     n124 = overlap_values$..count..[1] + overlap_values$..count..[5],
                     n134 = overlap_values$..count..[1] + overlap_values$..count..[3], 
                     n234 = overlap_values$..count..[1] + overlap_values$..count..[2], 
                     n1234 = overlap_values$..count..[1],
                     category = c("Females/LPS", "Females/VEH", "Males/LPS", "Males/VEH"), 
                     fill = c("red", "red", "skyblue1", "skyblue1"),
                     alpha = c(0.9, 0.2, 0.9, 0.4), cex = rep(2.5, 15), cat.cex = rep(2.5, 4),
                     cat.just = list(c(0.4,0), c(0.6,0), c(0.4,0), c(0.6,0))))
dev.off()

# 2) Test for enrichment of male eQTLs on x-chromosome ====

#Make 1x4 vectors of eQTL counts and merge them into a table
auto_vector <- c(nrow(auto_eQTLs[which(auto_eQTLs$MaleNull < 0.05),]),
                 nrow(auto_eQTLs[which(auto_eQTLs$MaleLPS < 0.05),]),
                 nrow(auto_eQTLs[which(auto_eQTLs$FemaleNull < 0.05),]),
                 nrow(auto_eQTLs[which(auto_eQTLs$FemaleLPS < 0.05),]))
x_vector <- c(nrow(x_eQTLs[which(x_eQTLs$MaleNull < 0.05),]),
                 nrow(x_eQTLs[which(x_eQTLs$MaleLPS < 0.05),]),
                 nrow(x_eQTLs[which(x_eQTLs$FemaleNull < 0.05),]),
                 nrow(x_eQTLs[which(x_eQTLs$FemaleLPS < 0.05),]))
test_table <- cbind.data.frame(auto_eQTLs=auto_vector, x_eQTLs=x_vector)
rownames(test_table) <- colnames(auto_eQTLs)
#Run tests
male_enrichment_4_way <- fisher.test(test_table)

#Make 1x2 vectors of eQTL counts and merge them into a table
auto_vector <- c(nrow(auto_eQTLs[which(auto_eQTLs$MaleNull < 0.05 | auto_eQTLs$MaleLPS < 0.05),]),
                 nrow(auto_eQTLs[which(auto_eQTLs$FemaleNull < 0.05 | auto_eQTLs$FemaleLPS < 0.05),]))
x_vector <- c(nrow(x_eQTLs[which(x_eQTLs$MaleNull < 0.05 | x_eQTLs$MaleLPS < 0.05),]),
                 nrow(x_eQTLs[which(x_eQTLs$FemaleNull < 0.05 | x_eQTLs$FemaleLPS < 0.05),]))
test_table <- cbind.data.frame(auto_eQTLs=auto_vector, x_eQTLs=x_vector)
rownames(test_table) <- c("Males", "Females")
#Run tests
male_enrichment_2_way <- fisher.test(test_table)

#Remove junk
remove(auto_vector, x_vector, test_table)

# 3) Create Export Lists for Enrichment Testing ====

#Make Lists
treatment_specific <-  mashr_raw[which(mashr_raw$MaleNull > 0.05 & mashr_raw$FemaleNull > 0.05 & 
                                         mashr_raw$MaleLPS < 0.05 & mashr_raw$FemaleLPS < 0.05), "gene"]
vehicle_specific <- mashr_raw[which(mashr_raw$MaleNull < 0.05 & mashr_raw$FemaleNull < 0.05 & 
                                      mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS > 0.05), "gene"]
male_specific <-  mashr_raw[which(mashr_raw$MaleNull < 0.05 & mashr_raw$FemaleNull > 0.05 & 
                                         mashr_raw$MaleLPS < 0.05 & mashr_raw$FemaleLPS > 0.05), "gene"]
female_specific <- mashr_raw[which(mashr_raw$MaleNull > 0.05 & mashr_raw$FemaleNull < 0.05 & 
                                      mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS < 0.05), "gene"]
female_LPS_specific <- mashr_raw[which(mashr_raw$MaleNull > 0.05 & mashr_raw$FemaleNull > 0.05 & 
                                         mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS < 0.05), "gene"]
male_LPS_specific <- mashr_raw[which(mashr_raw$MaleNull > 0.05 & mashr_raw$FemaleNull > 0.05 & 
                                         mashr_raw$MaleLPS < 0.05 & mashr_raw$FemaleLPS > 0.05), "gene"]
female_LPS_all <- mashr_raw[which(mashr_raw$FemaleLPS < 0.05), "gene"]
male_LPS_all <- mashr_raw[which(mashr_raw$MaleLPS < 0.05), "gene"]
female_VEH_all <- mashr_raw[which(mashr_raw$FemaleNull < 0.05), "gene"]
male_VEH_all <- mashr_raw[which(mashr_raw$MaleNull < 0.05), "gene"]

#Calculate numbers for manuscript
num_two_case <- length(treatment_specific) + length(vehicle_specific) + length(female_specific) + length(male_specific)
num_one_case <- length(female_LPS_specific) + length(male_LPS_specific) + length(mashr_raw[which(mashr_raw$MaleNull > 0.05 & mashr_raw$FemaleNull < 0.05 & mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS > 0.05), "gene"]) + length(mashr_raw[which(mashr_raw$MaleNull < 0.05 & mashr_raw$FemaleNull > 0.05 & mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS > 0.05), "gene"])

#Write gene names to files
write.table(treatment_specific, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/treatment_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(vehicle_specific, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/vehicle_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(male_specific, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/male_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(female_specific, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/female_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(female_LPS_specific, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/female_LPS_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(male_LPS_specific, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/male_LPS_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(female_LPS_all, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/female_LPS_all_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(male_LPS_all, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/male_LPS_all_genes.tsv", row.names = FALSE, quote = FALSE, col.names = FALSE)

# 4) Analyze Pre-Mashr eQTL Results ====

#Remove junk
suppressWarnings(remove(num_eQTLs, num_gemma_x_eQTLs, male_specific, female_specific, female_LPS_specific, male_LPS_specific, treatment_specific, vehicle_specific))

#Read in raw results
raw_gemma_eQTLs <- read.delim(file = "C:/Users/Janel/Documents/sex-specific/eQTL_analysis/combined.auto.gemma.assoc.tsv", 
                        stringsAsFactors = FALSE, sep = "\t")
#Filter out junk genes
raw_gemma_eQTLs <- raw_gemma_eQTLs[which(raw_gemma_eQTLs$gene %in% junk_genes == FALSE),]

#Filter for strongest eQTL for each gene
raw_gemma_eQTLs <- raw_gemma_eQTLs[order(raw_gemma_eQTLs$Min_p_score),]
gemma_eQTLs <- raw_gemma_eQTLs[duplicated(raw_gemma_eQTLs$gene) == FALSE,]
remove(raw_gemma_eQTLs)

#Adjust p-values
gemma_eQTLs$ML_p_score <- p.adjust(gemma_eQTLs$ML_p_score)
gemma_eQTLs$MV_p_score <- p.adjust(gemma_eQTLs$MV_p_score)
gemma_eQTLs$FL_p_score <- p.adjust(gemma_eQTLs$FL_p_score)
gemma_eQTLs$FV_p_score <- p.adjust(gemma_eQTLs$FV_p_score)

#Filter auto and x results 
gemma_auto_eQTLs <- gemma_eQTLs[which(gemma_eQTLs$gene %in% auto_genes),c("ML_p_score", "MV_p_score", "FL_p_score", "FV_p_score")]
gemma_x_eQTLs <- gemma_eQTLs[which(gemma_eQTLs$gene %in% x_genes),c("ML_p_score", "MV_p_score", "FL_p_score", "FV_p_score")]

#Rename columns
colnames(gemma_auto_eQTLs) <- c("MaleLPS", "MaleNull", "FemaleLPS", "FemaleNull")
colnames(gemma_x_eQTLs) <- c("MaleLPS", "MaleNull", "FemaleLPS", "FemaleNull")



#Count number of significant gene-snp combos in at least one condition
num_gemma_auto_eQTLs <- nrow(gemma_auto_eQTLs[which(gemma_auto_eQTLs$MaleNull < 0.05 | gemma_auto_eQTLs$MaleLPS < 0.05 | gemma_auto_eQTLs$FemaleNull < 0.05 | gemma_auto_eQTLs$FemaleLPS < 0.05),])
num_gemma_x_eQTLs <- nrow(gemma_x_eQTLs[which(gemma_x_eQTLs$MaleNull < 0.05 | gemma_x_eQTLs$MaleLPS < 0.05 | gemma_x_eQTLs$FemaleNull < 0.05 | gemma_x_eQTLs$FemaleLPS < 0.05),])
num_eQTLs <- num_gemma_auto_eQTLs + num_gemma_x_eQTLs

#Plot Autosomal eQTLs
#Calculate overlap groups
overlap_values<-get.venn.partitions(list("Females/LPS" = rownames(gemma_auto_eQTLs[which(gemma_auto_eQTLs$FemaleLPS < 0.05),]), 
                                         "Females/VEH" = rownames(gemma_auto_eQTLs[which(gemma_auto_eQTLs$FemaleNull < 0.05),]),
                                         "Males/LPS" = rownames(gemma_auto_eQTLs[which(gemma_auto_eQTLs$MaleLPS < 0.05),]), 
                                         "Males/VEH" = rownames(gemma_auto_eQTLs[which(gemma_auto_eQTLs$MaleNull < 0.05),])))
#Redefine output diagram for venn diagrams
out_dir <- "C:/Users/Janel/Documents/Sex-specific/Visuals/venn_diagrams/"
#Draw Venn Diagram
tiff(paste(out_dir, "gemma_auto_venn_diagram.jpg", sep = ""), width = 11.5, height = 8, units = 'in', res = 500)
print(draw.quad.venn(area1 = length(gemma_auto_eQTLs$FemaleLPS[which(gemma_auto_eQTLs$FemaleLPS < 0.05)]), 
                     area2 = length(gemma_auto_eQTLs$FemaleNull[which(gemma_auto_eQTLs$FemaleNull < 0.05)]), 
                     area3 = length(gemma_auto_eQTLs$MaleLPS[which(gemma_auto_eQTLs$MaleLPS < 0.05)]), 
                     area4 = length(gemma_auto_eQTLs$MaleNull[which(gemma_auto_eQTLs$MaleNull < 0.05)]), 
                     n12 = overlap_values$..count..[1] + overlap_values$..count..[5] + overlap_values$..count..[9] + overlap_values$..count..[13], 
                     n13 = overlap_values$..count..[1] + overlap_values$..count..[3] + overlap_values$..count..[9] + overlap_values$..count..[11], 
                     n14 = overlap_values$..count..[1] +  overlap_values$..count..[3] + overlap_values$..count..[5] + overlap_values$..count..[7], 
                     n23 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[9] + overlap_values$..count..[10], 
                     n24 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[5] + overlap_values$..count..[6], 
                     n34 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[3] + overlap_values$..count..[4], 
                     n123 = overlap_values$..count..[1] + overlap_values$..count..[9], 
                     n124 = overlap_values$..count..[1] + overlap_values$..count..[5],
                     n134 = overlap_values$..count..[1] + overlap_values$..count..[3], 
                     n234 = overlap_values$..count..[1] + overlap_values$..count..[2], 
                     n1234 = overlap_values$..count..[1],
                     category = c("Females/LPS", "Females/VEH", "Males/LPS", "Males/VEH"), 
                     fill = c("red", "red", "skyblue1", "skyblue1"),
                     alpha = c(0.9, 0.2, 0.9, 0.4), cex = rep(2.5, 15), cat.cex = rep(2.5, 4),
                     cat.just = list(c(0.4,0), c(0.6,0), c(0.4,0), c(0.6,0))))
dev.off()

# 5) Investigate Neutrophil Ontologies (Commented out in favor of iPathway Results)====

#Remove junk
suppressWarnings(remove(gemma_x_genes, gemma_auto_eQTLs, gemma_eQTLs, num_gemma_x_eQTLs, num_eQTLs, num_x_eQTLs, overlap_values, num_gemma_auto_eQTLs, x_genes, auto_genes))

#Read in table with neutrophil enrichment gene names
#female_go_enrichment <- read.delim(file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/Female_GO_enrichment.txt", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
#Filter for significant results
#female_go_enrichment <- female_go_enrichment[which(female_go_enrichment$Adjusted.P.value <= 0.05),]
#Get neutrophil rows and then gene list
#female_go_enrichment <- female_go_enrichment[grep("neutrophil", female_go_enrichment$Term),]
#female_neutrophil_genes <- intersect(str_split(female_go_enrichment[1,"Genes"], ";")[[1]], str_split(female_go_enrichment[2,"Genes"], ";")[[1]])
#female_neutrophil_genes <- intersect(female_neutrophil_genes, str_split(female_go_enrichment[3,"Genes"], ";")[[1]])
#Check for overlap with male_LPS_all gene list
#gene_overlap <- intersect(male_LPS_all, female_neutrophil_genes)
#non_overlap <- female_neutrophil_genes[which(female_neutrophil_genes %in% male_LPS_all == FALSE)]
#femaleLPS_exclusive <- non_overlap[which(non_overlap %in% mashr_raw[which(mashr_raw$MaleNull > 0.05 & mashr_raw$FemaleNull > 0.05 & 
#                                            mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS < 0.05), "gene"])]
#female_exclusive <- non_overlap[which(non_overlap %in% mashr_raw[which(mashr_raw$MaleNull > 0.05 & mashr_raw$FemaleNull < 0.05 & 
#                                                          mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS < 0.05), "gene"])]
#non_malesLPS <- non_overlap[which(non_overlap %in% mashr_raw[which(mashr_raw$MaleNull < 0.05 & mashr_raw$FemaleNull < 0.05 & 
#                                                                         mashr_raw$MaleLPS > 0.05 & mashr_raw$FemaleLPS < 0.05), "gene"])]
#non_malesLPS_20 <- intersect(female_neutrophil_genes, mashr_raw[which(mashr_raw$MaleLPS > 0.20), "gene"]) 
#femaleLPS_exclusive_20 <- intersect(female_neutrophil_genes, mashr_raw[which(mashr_raw$MaleLPS > 0.20 & mashr_raw$MaleNull > 0.20 & 
#                                                                               mashr_raw$FemaleNull > 0.20), "gene"])
#female_exclusive_20 <- intersect(female_neutrophil_genes, mashr_raw[which(mashr_raw$MaleLPS > 0.20 & mashr_raw$MaleNull > 0.20), "gene"])

# 6) Create files for iPathways analysis ====

#Remove junk
suppressWarnings(remove(female_exclusive, gene_overlap, female_neutrophil_genes, non_malesLPS, femaleLPS_exclusive, gemma_x_eQTLs))
suppressWarnings(remove(mashr_raw, x_eQTLs, female_go_enrichment, out_dir, non_malesLPS_20, female_exclusive_20, femaleLPS_exclusive_20))

#Read in Mashr output
Mash_Results_allgenes_strong_020320 <- readRDS("C:/Users/Janel/Documents/sex-specific/eQTL_analysis/Mash_Results_allgenes_strong_020320.rds")
#Pull posterior means and lfsrs
posterior_means <- as.data.frame(Mash_Results_allgenes_strong_020320$result$PosteriorMean)
posterior_SDs <- as.data.frame(Mash_Results_allgenes_strong_020320$result$PosteriorSD)
lfsrs <- as.data.frame(Mash_Results_allgenes_strong_020320$result$lfsr)
#Replace periods with colons
rownames(lfsrs) <- str_replace(rownames(lfsrs), "\\:", ":")
rownames(posterior_SDs) <- str_replace(rownames(posterior_SDs), "\\:.", ":")
rownames(posterior_means) <- str_replace(rownames(posterior_means), "\\:.", ":")
#Remove positions from rownames
rownames(lfsrs) <- substr(rownames(lfsrs),1,regexpr(":", rownames(lfsrs)) - 1)
rownames(posterior_SDs) <- substr(rownames(posterior_SDs),1,regexpr(":", rownames(posterior_SDs)) - 1)
rownames(posterior_means) <- substr(rownames(posterior_means),1,regexpr(":", rownames(posterior_means)) - 1)
#Remove junk genes
lfsrs <- lfsrs[which(rownames(lfsrs) %in% junk_genes == FALSE),]
posterior_SDs <- posterior_SDs[which(rownames(posterior_SDs) %in% junk_genes == FALSE),]
posterior_means <- posterior_means[which(rownames(posterior_means) %in% junk_genes == FALSE),]
#Create filters for lfsrs and posterior means and bind columns to make output files
femaleLPS <- cbind.data.frame(lfsr=lfsrs[female_LPS_all,"FemaleLPS"], posterior_mean=posterior_means[female_LPS_all,"FemaleLPS"])
rownames(femaleLPS) <- female_LPS_all
maleLPS <- cbind.data.frame(lfsr=lfsrs[male_LPS_all,"MaleLPS"], posterior_mean=posterior_means[male_LPS_all,"MaleLPS"])
rownames(maleLPS) <- male_LPS_all

#Write files for iPathways
write.table(femaleLPS, file = "C:/Users/Janel/Documents/sex-specific/eQTL_analysis/enrichment_results/female_LPS_ipath.tsv", row.names = TRUE, quote = FALSE, col.names = TRUE)
write.table(maleLPS, file = "C:/Users/Janel/Documents/sex-specific/eQTL_analysis/enrichment_results/male_LPS_ipath.tsv", row.names = TRUE, quote = FALSE, col.names = TRUE)
write.table(lfsrs, file = "C:/Users/Janel/Documents/sex-specific/eQTL_analysis/enrichment_results/lfsrs.tsv", row.names = TRUE, quote = FALSE, col.names = TRUE)
write.table(posterior_means, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/posterior_means.tsv", row.names = TRUE, quote = FALSE, col.names = TRUE)
write.table(posterior_SDs, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/enrichment_results/posterior_SDs.tsv", row.names = TRUE, quote = FALSE, col.names = TRUE)


# 7) Check for gene overlap with sex-by-treatment interaction list ====

#Read in sex-by-treatment interaction results
sex_by_treatment <- read.table(file = "C:/Users/Janel/Documents/sex-specific/overall_analyses/sbt.all.gemma.assoc.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE, row.names = 1)
colnames(sex_by_treatment) <- c("chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "p_score")
#Remove junk_genes
sex_by_treatment <- sex_by_treatment[which(rownames(sex_by_treatment) %in% junk_genes == FALSE),]
#Append p_adj_score
sex_by_treatment <- cbind.data.frame(sex_by_treatment, p_adj_score=p.adjust(sex_by_treatment$p_score, method = "fdr"))
#filter by pvals
sex_by_treatment <- sex_by_treatment[which(sex_by_treatment$p_adj_score < 0.05),]


#Get names of genes with eQTLs
eQTL_genes <- rownames(lfsrs[which(lfsrs$MaleNull < 0.05 | lfsrs$MaleLPS < 0.05 | lfsrs$FemaleNull < 0.05 | lfsrs$FemaleLPS < 0.05),])
#Calculate number of overlap genes
overlap <- intersect(eQTL_genes, rownames(sex_by_treatment))
overlap_four_conditions <- intersect(rownames(lfsrs[which(lfsrs$MaleNull < 0.05 & lfsrs$MaleLPS < 0.05 & lfsrs$FemaleNull < 0.05 & lfsrs$FemaleLPS < 0.05),]), rownames(sex_by_treatment))

#write overlap_four to file
write.csv(overlap_four_conditions, file = "C:/Users/Janel/Documents/sex-specific/eQTL_analysis/enrichment_results/overlap_four.csv", quote = FALSE, row.names = FALSE)

#Check directionality in various groups of overlap_four conditions list
overlap_four_conditions_directions <- as.data.frame(Mash_Results_allgenes_strong_020320$result$PosteriorMean)
rownames(overlap_four_conditions_directions) <- substr(rownames(overlap_four_conditions_directions),1,regexpr("\\:", rownames(overlap_four_conditions_directions)) - 1)
#Filter out junk genes
overlap_four_conditions_directions <- overlap_four_conditions_directions[which(rownames(overlap_four_conditions_directions) %in% junk_genes == FALSE),]
#Create list of genes with shared directionality in all four conditions irrespective of sbt status
shared_four_directions <- rownames(overlap_four_conditions_directions[which(rownames(overlap_four_conditions_directions) %in% eQTL_genes & 
                                                                              abs(abs(overlap_four_conditions_directions$MaleNull + 
                                                                                    overlap_four_conditions_directions$MaleLPS + 
                                                                                    overlap_four_conditions_directions$FemaleNull + 
                                                                                    overlap_four_conditions_directions$FemaleLPS) - 
                                                                              abs(overlap_four_conditions_directions$MaleNull) - 
                                                                              abs(overlap_four_conditions_directions$MaleLPS) - 
                                                                              abs(overlap_four_conditions_directions$FemaleNull) - 
                                                                              abs(overlap_four_conditions_directions$FemaleLPS)) < 0.00001),])
shared_four_conditions <- rownames(lfsrs[which(lfsrs$MaleNull < 0.05 & lfsrs$MaleLPS < 0.05 & lfsrs$FemaleNull < 0.05 & lfsrs$FemaleLPS < 0.05),])
shared_four_conditions_directions <- intersect(shared_four_conditions, shared_four_directions)

#Break up this block with a single line that makes matrix for overlap genes in addition to 4-condition overlap genes
overlap_directions <- overlap_four_conditions_directions[which(rownames(overlap_four_conditions_directions) %in% overlap),]
overlap_four_conditions_directions <- overlap_four_conditions_directions[which(rownames(overlap_four_conditions_directions) %in% overlap_four_conditions),]
overlap_four_conditions_directions <- cbind.data.frame(overlap_four_conditions_directions, all_four=
                                                         ifelse(abs(abs(overlap_four_conditions_directions$MaleNull + 
                                                                          overlap_four_conditions_directions$MaleLPS + 
                                                                          overlap_four_conditions_directions$FemaleNull + 
                                                                          overlap_four_conditions_directions$FemaleLPS) - 
                                                                      abs(overlap_four_conditions_directions$MaleNull) - 
                                                                      abs(overlap_four_conditions_directions$MaleLPS) - 
                                                                      abs(overlap_four_conditions_directions$FemaleNull) - 
                                                                      abs(overlap_four_conditions_directions$FemaleLPS)) < 0.00001, TRUE, FALSE))
#Check for enrichment of sbt results among 

#Check bidirectionality of sexes on overlap directions
overlap_directions <- cbind.data.frame(overlap_directions, opp_sexes=
                                                         ifelse(abs(overlap_directions$MaleNull + 
                                                                      overlap_directions$MaleLPS) == 
                                                                  abs(overlap_directions$MaleNull) + 
                                                                  abs(overlap_directions$MaleLPS) &
                                                                  abs(overlap_directions$FemaleNull + 
                                                                  overlap_directions$FemaleLPS) == 
                                                                  abs(overlap_directions$FemaleNull) + 
                                                                  abs(overlap_directions$FemaleLPS) & 
                                                                  abs(overlap_directions$MaleNull + 
                                                                        overlap_directions$MaleLPS + 
                                                                        overlap_directions$FemaleNull + 
                                                                        overlap_directions$FemaleLPS) != 
                                                                  abs(overlap_directions$MaleNull) + 
                                                                  abs(overlap_directions$MaleLPS) + 
                                                                  abs(overlap_directions$FemaleNull) + 
                                                                  abs(overlap_directions$FemaleLPS), TRUE, FALSE), 
                                       all_four=
                                         ifelse(abs(overlap_directions$MaleNull + 
                                                      overlap_directions$MaleLPS + 
                                                      overlap_directions$FemaleNull + 
                                                      overlap_directions$FemaleLPS) < 
                                                  abs(overlap_directions$MaleNull) + 
                                                  abs(overlap_directions$MaleLPS) + 
                                                  abs(overlap_directions$FemaleNull) + 
                                                  abs(overlap_directions$FemaleLPS), FALSE, TRUE))
#Make all directions
all_directions <- as.data.frame(Mash_Results_allgenes_strong_020320$result$PosteriorMean)
rownames(all_directions) <- substr(rownames(all_directions),1,regexpr("\\:", rownames(all_directions)) - 1)
all_directions <- all_directions[which(rownames(all_directions) %in% eQTL_genes & rownames(all_directions) %in% junk_genes == FALSE),]
#Check bidirectionality of sexes on all directions
all_directions <- cbind.data.frame(all_directions, opp_sexes=
                                         ifelse(abs(all_directions$MaleNull + 
                                                      all_directions$MaleLPS) == 
                                                  abs(all_directions$MaleNull) + 
                                                  abs(all_directions$MaleLPS) &
                                                  abs(all_directions$FemaleNull + 
                                                        all_directions$FemaleLPS) == 
                                                  abs(all_directions$FemaleNull) + 
                                                  abs(all_directions$FemaleLPS) & 
                                                  abs(all_directions$MaleNull + 
                                                        all_directions$MaleLPS + 
                                                        all_directions$FemaleNull + 
                                                        all_directions$FemaleLPS) != 
                                                  abs(all_directions$MaleNull) + 
                                                  abs(all_directions$MaleLPS) + 
                                                  abs(all_directions$FemaleNull) + 
                                                  abs(all_directions$FemaleLPS), TRUE, FALSE), 
                                                all_four=
                                                  ifelse(abs(all_directions$MaleNull + 
                                                               all_directions$MaleLPS + 
                                                               all_directions$FemaleNull + 
                                                               all_directions$FemaleLPS) < 
                                                           abs(all_directions$MaleNull) + 
                                                           abs(all_directions$MaleLPS) + 
                                                           abs(all_directions$FemaleNull) + 
                                                           abs(all_directions$FemaleLPS), FALSE, TRUE))

#Get significance for overlap
overlap_significance <- as.data.frame(Mash_Results_allgenes_strong_020320$result$lfsr)
rownames(overlap_significance) <- substr(rownames(overlap_significance),1,regexpr("\\:", rownames(overlap_significance)) - 1)
overlap_significance <- overlap_significance[which(rownames(overlap_significance) %in% overlap),]
#Remove junk genes
overlap_significance <- overlap_significance[which(rownames(overlap_significance) %in% junk_genes == FALSE),]
#Get all significance
all_significance <- lfsrs[rownames(all_directions),]
#Append columns to significance dfs
overlap_significance <- cbind.data.frame(overlap_significance, 
                                         all_four=ifelse(ifelse(overlap_significance$MaleNull<0.05, TRUE,FALSE) & 
                                                                                  ifelse(overlap_significance$MaleLPS<0.05, TRUE,FALSE) & 
                                                                                  ifelse(overlap_significance$FemaleNull<0.05, TRUE,FALSE) & 
                                                                                  ifelse(overlap_significance$FemaleLPS<0.05, TRUE,FALSE), TRUE, FALSE),
                                         one_sex=ifelse(ifelse(overlap_significance$MaleNull<0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$MaleLPS<0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleNull>0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleLPS>0.05, TRUE,FALSE) | 
                                                          ifelse(overlap_significance$MaleNull>0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$MaleLPS>0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleNull<0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleLPS<0.05, TRUE,FALSE), TRUE, FALSE),
                                         one_treatment=ifelse(ifelse(overlap_significance$MaleNull<0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleNull<0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$MaleLPS>0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleLPS>0.05, TRUE,FALSE) | 
                                                          ifelse(overlap_significance$MaleNull>0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleNull>0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$MaleLPS<0.05, TRUE,FALSE) & 
                                                          ifelse(overlap_significance$FemaleLPS<0.05, TRUE,FALSE), TRUE, FALSE))
all_significance <- cbind.data.frame(all_significance, 
                                         all_four=ifelse(ifelse(all_significance$MaleNull<0.05, TRUE,FALSE) & 
                                                           ifelse(all_significance$MaleLPS<0.05, TRUE,FALSE) & 
                                                           ifelse(all_significance$FemaleNull<0.05, TRUE,FALSE) & 
                                                           ifelse(all_significance$FemaleLPS<0.05, TRUE,FALSE), TRUE, FALSE),
                                         one_sex=ifelse(ifelse(all_significance$MaleNull<0.05, TRUE,FALSE) & 
                                                          ifelse(all_significance$MaleLPS<0.05, TRUE,FALSE) & 
                                                          ifelse(all_significance$FemaleNull>0.05, TRUE,FALSE) & 
                                                          ifelse(all_significance$FemaleLPS>0.05, TRUE,FALSE) | 
                                                          ifelse(all_significance$MaleNull>0.05, TRUE,FALSE) & 
                                                          ifelse(all_significance$MaleLPS>0.05, TRUE,FALSE) & 
                                                          ifelse(all_significance$FemaleNull<0.05, TRUE,FALSE) & 
                                                          ifelse(all_significance$FemaleLPS<0.05, TRUE,FALSE), TRUE, FALSE),
                                     one_treatment=ifelse(ifelse(all_significance$MaleNull<0.05, TRUE,FALSE) & 
                                                            ifelse(all_significance$FemaleNull<0.05, TRUE,FALSE) & 
                                                            ifelse(all_significance$MaleLPS>0.05, TRUE,FALSE) & 
                                                            ifelse(all_significance$FemaleLPS>0.05, TRUE,FALSE) | 
                                                            ifelse(all_significance$MaleNull>0.05, TRUE,FALSE) & 
                                                            ifelse(all_significance$FemaleNull>0.05, TRUE,FALSE) & 
                                                            ifelse(all_significance$MaleLPS<0.05, TRUE,FALSE) & 
                                                            ifelse(all_significance$FemaleLPS<0.05, TRUE,FALSE), TRUE, FALSE))

#Check percentages of four-way directional eGenes in overlap and all eGenes
nrow(overlap_directions[which(overlap_directions$all_four == TRUE),])/nrow(overlap_directions) #71.3%
nrow(all_directions[which(all_directions$all_four == TRUE),])/nrow(all_directions) #67.1%
fisher.test(matrix(data = c(nrow(overlap_directions[which(overlap_directions$all_four == TRUE),]),
                            nrow(overlap_directions)-nrow(overlap_directions[which(overlap_directions$all_four == TRUE),]),
                            nrow(all_directions[which(all_directions$all_four == TRUE),]),
                            nrow(all_directions)-nrow(all_directions[which(all_directions$all_four == TRUE),]))
                   , byrow = FALSE, nrow = 2, ncol = 2)) #p-value=0.2711
#Check for percentages of opposite directions of eGenes by sex in overlap and all eGenes
nrow(overlap_directions[which(overlap_directions$opp_sexes == TRUE),])/nrow(overlap_directions) #4.1%
nrow(all_directions[which(all_directions$opp_sexes == TRUE),])/nrow(all_directions) #9.8%
fisher.test(matrix(data = c(nrow(overlap_directions[which(overlap_directions$opp_sexes == TRUE),]),
                            nrow(overlap_directions)-nrow(overlap_directions[which(overlap_directions$opp_sexes == TRUE),]),
                            nrow(all_directions[which(all_directions$opp_sexes == TRUE),]),
                            nrow(all_directions)-nrow(all_directions[which(all_directions$opp_sexes == TRUE),]))
                   , byrow = FALSE, nrow = 2, ncol = 2)) #p-value=0.0131
#Check for percentages of four-condition eGenes in overlap and all eGenes
nrow(overlap_significance[which(overlap_significance$all_four == TRUE),])/nrow(overlap_significance) #17.5%
nrow(all_significance[which(all_significance$all_four == TRUE),])/nrow(all_significance) #23.5%
fisher.test(matrix(data = c(nrow(overlap_significance[which(overlap_significance$all_four == TRUE),]),
                            nrow(overlap_significance)-nrow(overlap_significance[which(overlap_significance$all_four == TRUE),]),
                            nrow(all_significance[which(all_significance$all_four == TRUE),]),
                            nrow(all_significance)-nrow(all_significance[which(all_significance$all_four == TRUE),]))
                   , byrow = FALSE, nrow = 2, ncol = 2)) #p-value=0.0889
#Check for percentages of sex-specific eGenes in overlap and all eGenes
nrow(overlap_significance[which(overlap_significance$one_sex == TRUE),])/nrow(overlap_significance) #5.3%
nrow(all_significance[which(all_significance$one_sex == TRUE),])/nrow(all_significance) #3.4%
fisher.test(matrix(data = c(nrow(overlap_significance[which(overlap_significance$one_sex == TRUE),]),
                            nrow(overlap_significance)-nrow(overlap_significance[which(overlap_significance$one_sex == TRUE),]),
                            nrow(all_significance[which(all_significance$one_sex == TRUE),]),
                            nrow(all_significance)-nrow(all_significance[which(all_significance$one_sex == TRUE),]))
                   , byrow = FALSE, nrow = 2, ncol = 2)) #p-value=0.198
#Check for percentages of treatment-specific eGenes in overlap and all eGenes
nrow(overlap_significance[which(overlap_significance$one_treatment == TRUE),])/nrow(overlap_significance) #1.2%
nrow(all_significance[which(all_significance$one_treatment == TRUE),])/nrow(all_significance) #1.5%
fisher.test(matrix(data = c(nrow(overlap_significance[which(overlap_significance$one_treatment == TRUE),]),
                            nrow(overlap_significance)-nrow(overlap_significance[which(overlap_significance$one_treatment == TRUE),]),
                            nrow(all_significance[which(all_significance$one_treatment == TRUE),]),
                            nrow(all_significance)-nrow(all_significance[which(all_significance$one_treatment == TRUE),]))
                   , byrow = FALSE, nrow = 2, ncol = 2)) #p-value>0.9999


# 8) Check for HLA Enrichment ====

#Remove some junk
suppressWarnings(remove(num_one_case, auto_eQTLs, num_two_case, overlap_four_conditions, male_enrichment_2_way, male_enrichment_4_way, maleLPS))

#Check for HLA Genes
HLA_eGenes <- eQTL_genes[which(regexpr("HLA", eQTL_genes) > 0)]

#Read in expr file and get list of expressed autosomal genes
expressed_genes <- read.table("C:/Users/Janel/Documents/sex-specific/processed_data/BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", header = TRUE, stringsAsFactors = FALSE)
#Make DE_genes to get list of all autosomal genes
DE_genes <- rownames(expressed_genes)
expressed_genes <- rownames(expressed_genes)
#Filter out junk genes
expressed_genes <- expressed_genes[which(expressed_genes %in% junk_genes == FALSE)]
#Check for actually being tested among eQTL genes
expressed_genes <- expressed_genes[which(expressed_genes %in% rownames(posterior_means))]
expressed_HLA <- expressed_genes[which(regexpr("HLA-", expressed_genes) > 0)]

#Check for enrichment among autosomal genes
HLA_table <- cbind.data.frame(HLA = c(length(HLA_eGenes), length(expressed_HLA) - length(HLA_eGenes)), ALL = c(num_auto_eQTLs, length(expressed_genes) - num_auto_eQTLs))
HLA_test <- fisher.test(HLA_table)

#HLA Venn Diagram
HLA_eQTLs <- lfsrs[expressed_HLA,]
#Append info about eGene status to HLA_eQTLs
HLA_expressed <- cbind.data.frame(HLA_eQTLs, eGene=ifelse(rownames(HLA_eQTLs) %in% HLA_eGenes, TRUE, FALSE))
HLA_eQTLs <- HLA_expressed[which(HLA_expressed$eGene == TRUE),]

#Calculate overlap groups
overlap_values<-get.venn.partitions(list("Females/LPS" = rownames(HLA_eQTLs[which(HLA_eQTLs$FemaleLPS < 0.05),]), 
                                         "Females/VEH" = rownames(HLA_eQTLs[which(HLA_eQTLs$FemaleNull < 0.05),]),
                                         "Males/LPS" = rownames(HLA_eQTLs[which(HLA_eQTLs$MaleLPS < 0.05),]), 
                                         "Males/VEH" = rownames(HLA_eQTLs[which(HLA_eQTLs$MaleNull < 0.05),])))
#Redefine output diagram for venn diagrams
out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/venn_diagrams/"
#Draw Venn Diagram
tiff(paste(out_dir, "mashr_HLA_venn_diagram.jpg", sep = ""), width = 11.5, height = 8, units = 'in', res = 500)
print(draw.quad.venn(area1 = length(HLA_eQTLs$FemaleLPS[which(HLA_eQTLs$FemaleLPS < 0.05)]), 
                     area2 = length(HLA_eQTLs$FemaleNull[which(HLA_eQTLs$FemaleNull < 0.05)]), 
                     area3 = length(HLA_eQTLs$MaleLPS[which(HLA_eQTLs$MaleLPS < 0.05)]), 
                     area4 = length(HLA_eQTLs$MaleNull[which(HLA_eQTLs$MaleNull < 0.05)]), 
                     n12 = overlap_values$..count..[1] + overlap_values$..count..[5] + overlap_values$..count..[9] + overlap_values$..count..[13], 
                     n13 = overlap_values$..count..[1] + overlap_values$..count..[3] + overlap_values$..count..[9] + overlap_values$..count..[11], 
                     n14 = overlap_values$..count..[1] +  overlap_values$..count..[3] + overlap_values$..count..[5] + overlap_values$..count..[7], 
                     n23 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[9] + overlap_values$..count..[10], 
                     n24 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[5] + overlap_values$..count..[6], 
                     n34 = overlap_values$..count..[1] + overlap_values$..count..[2] + overlap_values$..count..[3] + overlap_values$..count..[4], 
                     n123 = overlap_values$..count..[1] + overlap_values$..count..[9], 
                     n124 = overlap_values$..count..[1] + overlap_values$..count..[5],
                     n134 = overlap_values$..count..[1] + overlap_values$..count..[3], 
                     n234 = overlap_values$..count..[1] + overlap_values$..count..[2], 
                     n1234 = overlap_values$..count..[1],
                     category = c("Females/LPS", "Females/VEH", "Males/LPS", "Males/VEH"), 
                     fill = c("red", "red", "skyblue1", "skyblue1"),
                     alpha = c(0.9, 0.2, 0.9, 0.4), cex = rep(2.5, 15), cat.cex = rep(2.5, 4),
                     cat.just = list(c(0.4,0), c(0.6,0), c(0.4,0), c(0.6,0))))
dev.off()
#Clean-up
remove(out_dir, overlap_values)

#Check HLA directionalities 
HLA_directions <- posterior_means[HLA_eGenes,]
HLA_eQTLs <- cbind.data.frame(HLA_eQTLs, all_four_directions=
                                ifelse(abs(HLA_directions$MaleNull + 
                                             HLA_directions$MaleLPS + 
                                             HLA_directions$FemaleNull + 
                                             HLA_directions$FemaleLPS) < 
                                         abs(HLA_directions$MaleNull) + 
                                         abs(HLA_directions$MaleLPS) + 
                                         abs(HLA_directions$FemaleNull) + 
                                         abs(HLA_directions$FemaleLPS), FALSE, TRUE))
HLA_directions <- posterior_means[rownames(HLA_expressed),]
HLA_expressed <- cbind.data.frame(HLA_expressed, all_four_directions=
                                ifelse(abs(HLA_directions$MaleNull + 
                                             HLA_directions$MaleLPS + 
                                             HLA_directions$FemaleNull + 
                                             HLA_directions$FemaleLPS) < 
                                         abs(HLA_directions$MaleNull) + 
                                         abs(HLA_directions$MaleLPS) + 
                                         abs(HLA_directions$FemaleNull) + 
                                         abs(HLA_directions$FemaleLPS), FALSE, TRUE))
remove(HLA_directions)

#Check overlap with sex-by-treatment list
#Read in sex-by-treatment interaction results
sex_by_treatment <- read.table(file = "C:/Users/Janel/Documents/Sex-specific/overall_analyses/sbt.all.gemma.assoc.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(sex_by_treatment) <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "p_score")
#Remove junk genes from sex_by_treatment
sex_by_treatment <- sex_by_treatment[which(sex_by_treatment$gene %in% junk_genes == FALSE),]
#Append p_adj_score column
sex_by_treatment <- cbind.data.frame(sex_by_treatment, p_adj_score=p.adjust(sex_by_treatment$p_score, "fdr"))
#Sort data frames to match order
HLA_eQTLs <- HLA_eQTLs[order(rownames(HLA_eQTLs)),]
HLA_expressed <- HLA_expressed[order(rownames(HLA_expressed)),]
sex_by_treatment <- sex_by_treatment[order(sex_by_treatment$gene),]
#Append columns
HLA_eQTLs <- cbind.data.frame(HLA_eQTLs, 
                              sbt_adj_p_score=sex_by_treatment[which(sex_by_treatment$gene %in% rownames(HLA_eQTLs)),"p_adj_score"], 
                              sbt=ifelse(rownames(HLA_eQTLs) %in% sex_by_treatment[sex_by_treatment$p_adj_score <0.05,"gene"], TRUE, FALSE))
HLA_expressed <- cbind.data.frame(HLA_expressed, 
                              sbt_adj_p_score=sex_by_treatment[which(sex_by_treatment$gene %in% rownames(HLA_expressed)),"p_adj_score"], 
                              sbt=ifelse(rownames(HLA_expressed) %in% sex_by_treatment[sex_by_treatment$p_adj_score <0.05,"gene"], TRUE, FALSE))
#Count total sex-by-treatment genes
total_sbt <- nrow(sex_by_treatment[which(sex_by_treatment$p_adj_score < 0.05),])
#Clean-up
remove(sex_by_treatment)


#Check for eQTL enrichment of sbt results
sbt_HLA_eGene_table <- cbind.data.frame(HLA = c(length(HLA_eGenes) - nrow(HLA_eQTLs[which(HLA_eQTLs$sbt == TRUE),]), nrow(HLA_eQTLs[which(HLA_eQTLs$sbt == TRUE),])), 
                                  ALL = c(length(eQTL_genes[eQTL_genes %in% expressed_genes]) - length(overlap), length(overlap)))
sbt_HLA_eGene_test <- fisher.test(sbt_HLA_eGene_table)
#Check for enrichment of sbt results among all HLA genes relative to all 12,478 tested genes
sbt_HLA_expressed_table <- cbind.data.frame(HLA = c(length(expressed_HLA) - nrow(HLA_expressed[which(HLA_expressed$sbt == TRUE),]), nrow(HLA_expressed[which(HLA_expressed$sbt == TRUE),])), 
                                        ALL = c(length(DE_genes) - total_sbt, total_sbt))
sbt_HLA_expressed_test <- fisher.test(sbt_HLA_expressed_table)

#sort posterior means and lfsrs and filter the same way to prepare for tests
posterior_means <- posterior_means[order(rownames(posterior_means)),]
lfsrs <- lfsrs[order(rownames(lfsrs)),]

#Check for enrichment of 4-way shared results (all eGene denominators)
four_way_HLA_table_one <- cbind.data.frame(HLA = c(nrow(HLA_eQTLs[which(HLA_eQTLs$MaleNull < 0.05 & HLA_eQTLs$MaleLPS < 0.05 &
                                                                          HLA_eQTLs$FemaleNull < 0.05 & HLA_eQTLs$FemaleLPS < 0.05 &
                                                                          HLA_eQTLs$all_four_directions == TRUE),]), 
                                               length(HLA_eGenes) - nrow(HLA_eQTLs[which(HLA_eQTLs$MaleNull < 0.05 & HLA_eQTLs$MaleLPS < 0.05 &
                                                                                           HLA_eQTLs$FemaleNull < 0.05 & HLA_eQTLs$FemaleLPS < 0.05 &
                                                                                           HLA_eQTLs$all_four_directions == TRUE),])), 
                                    ALL = c(nrow(lfsrs[which(lfsrs$MaleNull < 0.05 & lfsrs$MaleLPS < 0.05 & 
                                                               lfsrs$FemaleNull < 0.05 & lfsrs$FemaleLPS < 0.05 &
                                                               rownames(lfsrs) %in% expressed_genes & 
                                                               abs(abs(posterior_means$MaleNull + 
                                                                     posterior_means$MaleLPS + 
                                                                     posterior_means$FemaleNull + 
                                                                     posterior_means$FemaleLPS) - 
                                                               (abs(posterior_means$MaleNull) + 
                                                               abs(posterior_means$MaleLPS) + 
                                                               abs(posterior_means$FemaleNull) + 
                                                               abs(posterior_means$FemaleLPS))) < 0.000001
                                                               ),]),
                                            nrow(lfsrs[which((lfsrs$MaleNull < 0.05 | lfsrs$MaleLPS < 0.05 | 
                                                               lfsrs$FemaleNull < 0.05 | lfsrs$FemaleLPS < 0.05) &
                                                               rownames(lfsrs) %in% expressed_genes),]) -
                                      nrow(lfsrs[which(lfsrs$MaleNull < 0.05 & lfsrs$MaleLPS < 0.05 & 
                                                         lfsrs$FemaleNull < 0.05 & lfsrs$FemaleLPS < 0.05 &
                                                         rownames(lfsrs) %in% expressed_genes & 
                                                         abs(abs(posterior_means$MaleNull + 
                                                               posterior_means$MaleLPS + 
                                                               posterior_means$FemaleNull + 
                                                               posterior_means$FemaleLPS) - 
                                                         (abs(posterior_means$MaleNull) + 
                                                         abs(posterior_means$MaleLPS) + 
                                                         abs(posterior_means$FemaleNull) + 
                                                         abs(posterior_means$FemaleLPS))) < 0.00001
                                      ),])))
four_way_HLA_test_one <- fisher.test(four_way_HLA_table_one)


#Check for enrichment of 4-way shared results (4-way share eGene denominators)
four_way_HLA_table_two <- cbind.data.frame(HLA = c(nrow(HLA_eQTLs[which(HLA_eQTLs$MaleNull < 0.05 & HLA_eQTLs$MaleLPS < 0.05 &
                                                                          HLA_eQTLs$FemaleNull < 0.05 & HLA_eQTLs$FemaleLPS < 0.05 &
                                                                          HLA_eQTLs$all_four_directions == TRUE),]), 
                                                   nrow(HLA_eQTLs[which(HLA_eQTLs$MaleNull < 0.05 & HLA_eQTLs$MaleLPS < 0.05 &
                                                                          HLA_eQTLs$FemaleNull < 0.05 & HLA_eQTLs$FemaleLPS < 0.05 ),]) - 
                                                     nrow(HLA_eQTLs[which(HLA_eQTLs$MaleNull < 0.05 & HLA_eQTLs$MaleLPS < 0.05 &
                                                                            HLA_eQTLs$FemaleNull < 0.05 & HLA_eQTLs$FemaleLPS < 0.05 &
                                                                            HLA_eQTLs$all_four_directions == TRUE),])), 
                                           ALL = c(nrow(lfsrs[which(lfsrs$MaleNull < 0.05 & lfsrs$MaleLPS < 0.05 & 
                                                                      lfsrs$FemaleNull < 0.05 & lfsrs$FemaleLPS < 0.05 &
                                                                      rownames(lfsrs) %in% expressed_genes & 
                                                                      abs(abs(posterior_means$MaleNull + 
                                                                            posterior_means$MaleLPS + 
                                                                            posterior_means$FemaleNull + 
                                                                            posterior_means$FemaleLPS) -  
                                                                      (abs(posterior_means$MaleNull) + 
                                                                      abs(posterior_means$MaleLPS) + 
                                                                      abs(posterior_means$FemaleNull) + 
                                                                      abs(posterior_means$FemaleLPS))) < 0.00001
                                           ),]),
                                           nrow(lfsrs[which(lfsrs$MaleNull < 0.05 & lfsrs$MaleLPS < 0.05 & 
                                                              lfsrs$FemaleNull < 0.05 & lfsrs$FemaleLPS < 0.05 &
                                                              rownames(lfsrs) %in% expressed_genes),]) - 
                                             nrow(lfsrs[which(lfsrs$MaleNull < 0.05 & lfsrs$MaleLPS < 0.05 & 
                                                                lfsrs$FemaleNull < 0.05 & lfsrs$FemaleLPS < 0.05 &
                                                                rownames(lfsrs) %in% expressed_genes & 
                                                                abs(abs(posterior_means$MaleNull + 
                                                                      posterior_means$MaleLPS + 
                                                                      posterior_means$FemaleNull + 
                                                                      posterior_means$FemaleLPS) - 
                                                                (abs(posterior_means$MaleNull) + 
                                                                abs(posterior_means$MaleLPS) + 
                                                                abs(posterior_means$FemaleNull) + 
                                                                abs(posterior_means$FemaleLPS))) < 0.00001
                                             ),])))
four_way_HLA_test_two <- fisher.test(four_way_HLA_table_two)


# 9) Output HLA table ====

#Remake HLA_eQTLs
HLA_output <- cbind.data.frame(lfsrs[HLA_eGenes,], posterior_means[HLA_eGenes,])
colnames(HLA_output) <- append(paste(colnames(lfsrs), "_lfsrs", sep = ""), paste(colnames(posterior_means), "_post_means", sep = ""))
#Sort column names
HLA_output <- HLA_output[,order(colnames(HLA_output))]
#Write to file
write.table(HLA_output, file = "C:/Users/Janel/Documents/Sex-specific/eQTL_analysis/supp_table_6_prep.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


