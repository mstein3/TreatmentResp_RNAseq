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
library(Hmisc)


# 1) Read in files for X-linked-GRM df change ====

#Set input directory
inp_dir <- "C:/Users/Janel/Documents/sex-specific/processed_data/"
filter_dir <- "C:/Users/Janel/Documents/sex-specific/filters/"

#Define genders
genders <- c("females", "males")

#Create empty data frames to be appended to (number of rows manually obtained by looking at file)
expr_final <- as.data.frame(matrix(nrow = 393))
covar <- as.data.frame(matrix(ncol = 7))
colnames(covar) <- c("Extract_Ord", "sex", "treat", "Age", "Pool_1.0", "library_tech", "FINDIV")
foldchange_df<-as.data.frame(matrix(nrow = 393))
foldchange_df_x<-as.data.frame(matrix(nrow = 393))

#Repeat for both genders
for (gender in genders) {
  
  #Read in expression and sort
  expr_raw <- read.table(paste(inp_dir, gender, ".logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  expr_raw <- expr_raw[order(rownames(expr_raw)),]
  
  ## Next block verifies genes are being added together correctly ##
  #Check whether on first or second pass of loop
  if (rownames(expr_final)[1] == "1") {
    rownames(expr_final) <- rownames(expr_raw)
  } else {
    expr_final <- expr_final[which(rownames(expr_final) %in% rownames(expr_raw)),]
    foldchange_df <- foldchange_df[which(rownames(foldchange_df) %in% rownames(expr_final)),]
    foldchange_df_x <- foldchange_df[which(rownames(foldchange_df) %in% rownames(expr_final)),]
  }
  #Add to expression final
  expr_raw <- expr_raw[which(rownames(expr_raw) %in% rownames(expr_final)),]
  expr_final <- cbind(expr_final, expr_raw)
  
  #Read in pheno file
  pheno_raw <- read.table(paste(inp_dir, gender, ".Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  #Remove junk from pheno file (need to read in id list from filter_dir)
  iid_list <- read.table(paste(filter_dir, gender, ".iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  pheno_raw <- pheno_raw[which(pheno_raw$FINDIV %in% iid_list$V2),]
  remove(iid_list)
  
  #Get list of ids
  individuals<-pheno_raw[,"FINDIV"]
  individuals<-subset(individuals,duplicated(individuals)!=TRUE)
  individuals_bound<-cbind(individuals,individuals,individuals)
  for (k in 1:nrow(individuals_bound)) {
    individuals_bound[k,2]<-pheno_raw[which(pheno_raw$FINDIV==individuals_bound[k,1] & pheno_raw$treat=="LPS"),][,"Extract_Ord"]
    individuals_bound[k,3]<-pheno_raw[which(pheno_raw$FINDIV==individuals_bound[k,1] & pheno_raw$treat=="null"),][,"Extract_Ord"]
  }
  remove(individuals,k)
  colnames(individuals_bound) <- c("FINDIV", "LPS", "null")
  
  #Filter individuals_bound for the X_GRM (need to read in X_GRM iids file)
  x_grm_ids <- read.table(paste(filter_dir, gender, ".X.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  individuals_bound_x <- individuals_bound[which(individuals_bound[,2] %in% x_grm_ids[,2]),]
  
  #Define fold change column
  foldchange<-vector()
  foldchange_x<-vector()
  #auto_grm
  for (i in 1:nrow(expr_raw)) {
    
    #Create fold change column (method 1 divide then average)
    templist<-vector()
    for (j in 1:nrow(individuals_bound)) {
      templist[j]<-expr_raw[i,sprintf("X%04d", individuals_bound[j,2])] - expr_raw[i,sprintf("X%04d", individuals_bound[j,3])]
    }
    #average values
    foldchange[i]<-mean(templist)
    
  }
  suppressWarnings(remove(i,j))
  #Confirm rownames of 
  #cbind foldchange to df
  foldchange_df<-cbind(foldchange_df, foldchange)
  rownames(foldchange_df) <- rownames(expr_raw)
  #x_grm
  for (i in 1:nrow(expr_raw)) {
    
    #Create fold change column (method 1 divide then average)
    templist<-vector()
    for (j in 1:nrow(individuals_bound_x)) {
      templist[j]<-expr_raw[i,sprintf("X%04d", individuals_bound_x[j,2])] - expr_raw[i,sprintf("X%04d", individuals_bound_x[j,3])]
    }
    #average values
    foldchange_x[i]<-mean(templist)
    
  }
  suppressWarnings(remove(i,j))
  #cbind foldchange to df
  foldchange_df_x<-cbind(foldchange_df_x, foldchange_x)
  rownames(foldchange_df_x) <- rownames(expr_raw)
  
  #Append pheno info to covar file
  covar <- rbind(covar, pheno_raw[, c("Extract_Ord", "sex", "treat", "Age", "Pool_1.0", "library_tech", "FINDIV")])
  
  #Sort data frames
  foldchange_df <- foldchange_df[order(row.names(foldchange_df)),]
  foldchange_df_x <- foldchange_df_x[order(row.names(foldchange_df_x)),]
  expr_final <- expr_final[order(row.names(expr_final)),]
}

#Remove trash from data frames
covar <- covar[2:nrow(covar),]
expr_final <- expr_final[,2:ncol(expr_final)]
foldchange_df <- foldchange_df[,2:ncol(foldchange_df)]
foldchange_df_x <- foldchange_df_x[,2:ncol(foldchange_df_x)]

#Define columns and rows of foldchange_df
colnames(foldchange_df)<-paste(genders,"_FC",sep = "")
rownames(foldchange_df)<-rownames(expr_raw)
colnames(foldchange_df_x)<-paste(genders,"_FC",sep = "")
rownames(foldchange_df_x)<-rownames(expr_raw)

#Clear junk
suppressWarnings(remove(pheno_raw, expr_raw, inp_dir,templist, foldchange, foldchange_x, x_grm_ids))

#Convert sex to words
covar[,"sex"] <- as.character(covar[,"sex"])
for (i in 1:nrow(covar)) {
  if (covar[i,2] == "0") {
    covar[i,2] <- "Males"
  } else if(covar[i,2] == "1"){
    covar[i,2] <- "Females" 
  }
}

#Clear junk
suppressWarnings(remove(covar, expr_final, i, gender, individuals_bound, individuals_bound_x))

# 2) Read in files for autosome-GRM df change ====

#Set input directory
inp_dir <- "C:/Users/Janel/Documents/sex-specific/processed_data/"
filter_dir <- "C:/Users/Janel/Documents/sex-specific/filters/"

#Define genders
genders <- c("females", "males")

#Create empty data frames to be appended to (Get row numbers from manual check of output file)
expr_final <- as.data.frame(matrix(nrow = 12478))
covar <- as.data.frame(matrix(ncol = 7))
colnames(covar) <- c("Extract_Ord", "sex", "treat", "Age", "Pool_1.0", "library_tech", "FINDIV")
foldchange_df_auto<-as.data.frame(matrix(nrow = 12478))

#Repeat for both genders
for (gender in genders) {
  
  #Read in id files
  iid_file <- read.table(paste(filter_dir, gender, ".iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  HLA_iid_file <- read.table(paste(filter_dir, "HLA.", gender, ".iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  #Read in list of HLA genes
  HLA_genes <- read.table(paste(filter_dir, "HLA_auto_genes_DE.tsv", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  HLA_genes <- as.character(HLA_genes$V1)
  
  #Read in expression
  expr_raw <- read.table(paste(inp_dir, "BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  #sort by row names
  expr_raw <- expr_raw[order(row.names(expr_raw)),]
  
  #reassign expr_raw to expr_final since this expression matrix has both genders
  expr_final<-expr_raw
  
  #Read in pheno file
  pheno_raw <- read.table(paste(inp_dir, gender, ".Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  #Filter pheno file for gender ids
  pheno_raw <- pheno_raw[which(pheno_raw$FINDIV %in% iid_file$V2),]
  remove(iid_file)
  #Filter pheno file for HLA ids
  pheno_HLA <- pheno_raw[which(pheno_raw$FINDIV %in% HLA_iid_file$V2),]
  remove(HLA_iid_file)
  
  #Get list of ids
  individuals<-pheno_raw[,"FINDIV"]
  individuals<-subset(individuals,duplicated(individuals)!=TRUE)
  individuals_bound<-cbind(individuals,individuals,individuals)
  for (k in 1:nrow(individuals_bound)) {
    individuals_bound[k,2]<-pheno_raw[which(pheno_raw$FINDIV==individuals_bound[k,1] & pheno_raw$treat=="LPS"),][,"Extract_Ord"]
    individuals_bound[k,3]<-pheno_raw[which(pheno_raw$FINDIV==individuals_bound[k,1] & pheno_raw$treat=="null"),][,"Extract_Ord"]
  }
  remove(individuals,k)
  colnames(individuals_bound) <- c("FINDIV", "LPS", "null")
  
  #Get list of HLA_ids
  HLA_individuals<-pheno_HLA[,"FINDIV"]
  HLA_individuals<-subset(HLA_individuals,duplicated(HLA_individuals)!=TRUE)
  HLA_individuals_bound<-cbind(HLA_individuals,HLA_individuals,HLA_individuals)
  for (k in 1:nrow(HLA_individuals_bound)) {
    HLA_individuals_bound[k,2]<-pheno_HLA[which(pheno_HLA$FINDIV==HLA_individuals_bound[k,1] & pheno_HLA$treat=="LPS"),][,"Extract_Ord"]
    HLA_individuals_bound[k,3]<-pheno_HLA[which(pheno_HLA$FINDIV==HLA_individuals_bound[k,1] & pheno_HLA$treat=="null"),][,"Extract_Ord"]
  }
  remove(HLA_individuals,k)
  colnames(HLA_individuals_bound) <- c("FINDIV", "LPS", "null")
  
  #Define fold change column
  foldchange_auto<-vector()
  #auto_grm
  for (i in 1:nrow(expr_raw)) {
    #Check whether gene is one of the special HLA_genes
    if(rownames(expr_raw)[i] %in% HLA_genes){
      #Create fold change column (method 1 divide then average)
      templist<-vector()
      for (j in 1:nrow(HLA_individuals_bound)) {
        templist[j]<-expr_raw[i,paste("X",HLA_individuals_bound[j,2],sep = "")] - expr_raw[i,paste("X",HLA_individuals_bound[j,3],sep = "")]
      }
      #average values
      foldchange_auto[i]<-mean(templist)
      
    } else {
      #Create fold change column (method 1 divide then average)
      templist<-vector()
      for (j in 1:nrow(individuals_bound)) {
        templist[j]<-expr_raw[i,paste("X",individuals_bound[j,2],sep = "")] - expr_raw[i,paste("X",individuals_bound[j,3],sep = "")]
      }
      #average values
      foldchange_auto[i]<-mean(templist)
    }
    
  }
  suppressWarnings(remove(i,j))
  #Confirm rownames of 
  #cbind foldchange to df
  foldchange_df_auto<-cbind(foldchange_df_auto, foldchange_auto)
  rownames(foldchange_df_auto) <- rownames(expr_raw)
  
  #Append pheno info to covar file
  covar <- rbind(covar, pheno_raw[, c("Extract_Ord", "sex", "treat", "Age", "Pool_1.0", "library_tech", "FINDIV")])
  
  #Sort data frames
  foldchange_df_auto <- foldchange_df_auto[order(row.names(foldchange_df_auto)),]
  expr_final <- expr_final[order(row.names(expr_final)),]
}

#Remove trash from data frames
covar <- covar[2:nrow(covar),]
foldchange_df_auto <- foldchange_df_auto[,2:ncol(foldchange_df_auto)]

#Define columns and rows of foldchange_df_auto
colnames(foldchange_df_auto)<-paste(genders,"_FC",sep = "")
rownames(foldchange_df_auto)<-rownames(expr_raw)

#Clear junk
suppressWarnings(remove(pheno_raw, expr_raw, inp_dir,templist, foldchange_auto, foldchange_x_auto, x_grm_ids, individuals_bound_x, individuals_bound))

#Convert sex to words
covar[,"sex"] <- as.character(covar[,"sex"])
for (i in 1:nrow(covar)) {
  if (covar[i,2] == "0") {
    covar[i,2] <- "Males"
  } else if(covar[i,2] == "1"){
    covar[i,2] <- "Females" 
  }
}

#Clear junk
remove(covar, expr_final, i)

#Remove junk genes
junk_genes <- c("HLA-DPA1", "HLA-DRB5", "HLA-G")
foldchange_df_auto <- foldchange_df_auto[rownames(foldchange_df_auto) %in% junk_genes == FALSE,]

# 3) Read in p-values ====

#Order of these next two lines is important!!!!!
#Merge x-GRM for x-linked genes with auto-GRM for autosomal genes into one data frame
foldchange_df_x <- rbind(cbind(foldchange_df_x, "Chr_type"=rep("X-linked", nrow(foldchange_df_x))), cbind(foldchange_df_auto, "Chr_type"=rep("Autosomal", nrow(foldchange_df_auto))))
#Merge two autosomal-GRM-based foldchange_dfs into one df
foldchange_df <- rbind(cbind(foldchange_df, "Chr_type"=rep("X-linked", nrow(foldchange_df))), cbind(foldchange_df_auto, "Chr_type"=rep("Autosomal", nrow(foldchange_df_auto))))


#Set input directory
inp_dir <- "C:/Users/Janel/Documents/sex-specific/de_analysis/"
process_dir <- "C:/Users/Janel/Documents/sex-specific/processed_data/"
out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/"

#make copy of foldchange df to hold significance data
scatter_df <- foldchange_df
scatter_df_x <- foldchange_df_x

#Loop through both genders
for (gender in genders) {
  
  ## Manipulate X-linked files##
  #read in x-linked files
  x_raw <- read.table(paste(inp_dir, gender, ".all.X.gemma.assoc.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  auto_raw <- read.table(paste(inp_dir, gender, ".all.auto.gemma.assoc.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  #Define colnames
  colnames(x_raw) <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score")
  colnames(auto_raw) <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score")
  #Convert final columns to numeric data
  x_raw[,"p_score"] <- as.numeric(x_raw[,"p_score"])
  auto_raw[,"p_score"] <- as.numeric(auto_raw[,"p_score"])
  x_raw[,"beta"] <- as.numeric(x_raw[,"beta"])
  auto_raw[,"beta"] <- as.numeric(auto_raw[,"beta"])
  #Sort data frames
  attach(x_raw)
  x_raw <- x_raw[order(gene),]
  detach(x_raw)
  attach(auto_raw)
  auto_raw <- auto_raw[order(gene),]
  detach(auto_raw)

  #Make rownames gene names for both data frames too
  rownames(auto_raw) <- auto_raw$gene
  rownames(x_raw) <- x_raw$gene
  
  #Remove junk genes
  auto_raw <- auto_raw[rownames(auto_raw) %in% junk_genes == FALSE,]
  
  #Merge auto_raw and auto_raw_auto into single dataframes for each GRM
  #i.e. merge x-linked and autosomal genes together
  auto_raw_x <- auto_raw
  auto_raw_x[rownames(x_raw),] <- x_raw
  
  #Remove junk
  suppressWarnings(remove(x_raw))
  
  #Merge Pvals into a df and adjust pvals
  #auto first
  pvals <- as.data.frame(cbind(auto_raw$gene, auto_raw$p_score, as.character(auto_raw$beta)))
  colnames(pvals) <- c("gene", "p_score", "beta")
  pvals[,"p_score"] <- as.numeric(p.adjust(auto_raw$p_score, method = "fdr"))
  pvals$auto_beta <- as.numeric(as.character(pvals$beta))
  #then x
  pvals_x <- as.data.frame(cbind(auto_raw_x$gene, auto_raw_x$p_score, as.character(auto_raw_x$beta)))
  colnames(pvals_x) <- c("gene", "p_score", "beta")
  pvals_x[,"p_score"] <- as.numeric(p.adjust(auto_raw_x$p_score, method = "fdr"))
  pvals_x$auto_beta <- as.numeric(as.character(pvals_x$beta))
  
  #Make column for significance
  #auto first
  pvals <- cbind(pvals, sig=ifelse(pvals[,2]<=0.05, "significant", "not significant"))
  #Convert to sig to factor
  pvals[,"sig"] <- as.factor(pvals[,"sig"])
  #then x
  pvals_x <- cbind(pvals_x, sig=ifelse(pvals_x[,2]<=0.05, "significant", "not significant"))
  #Convert to sig to factor
  pvals_x[,"sig"] <- as.factor(pvals_x[,"sig"])
  
  #Add -log columns
  pvals <- cbind(pvals, neg_log_auto=-log(pvals[,"p_score"], base = 10))
  pvals_x <- cbind(pvals_x, neg_log_auto=-log(pvals_x[,"p_score"], base = 10))
  
  #Order foldchange_df and pvals
  pvals <- pvals[order(pvals$gene),]
  pvals_x <- pvals_x[order(pvals_x$gene),]
  foldchange_df <- foldchange_df[order(row.names(foldchange_df)),]
  foldchange_df_x <- foldchange_df_x[order(row.names(foldchange_df_x)),]

  
  #Combine foldchanges and pvals into one df
  #auto first
  combined <- cbind(pvals[which(pvals$gene %in% rownames(foldchange_df)),], foldchange_df[which(rownames(foldchange_df) %in% pvals$gene),])
  colnames(combined) <- append(colnames(pvals), c("auto_females_FC", "auto_males_FC", "Chr_type"))
  #add significance columns to combined
  combined <- cbind(combined, auto_sig=ifelse(combined$p_score <= 0.05, "adj.P.Val<0.05", "Not Sig"))
  #Then x
  combined_x <- cbind(pvals_x[which(pvals_x$gene %in% rownames(foldchange_df_x)),], foldchange_df_x[which(rownames(foldchange_df_x) %in% pvals_x$gene),])
  colnames(combined_x) <- append(colnames(pvals_x), c("auto_females_FC", "auto_males_FC", "Chr_type"))
  #add significance columns to combined
  combined_x <- cbind(combined_x, auto_sig=ifelse(combined_x$p_score <= 0.05, "adj.P.Val<0.05", "Not Sig"))
  
  #filter scatter_df for rows in combined df and vice versa
  #auto first
  scatter_df <- scatter_df[order(rownames(scatter_df)),]
  scatter_df <- scatter_df[which(rownames(scatter_df) %in% combined$gene),]
  combined <- combined[which(combined$gene %in% rownames(scatter_df)),]
  #x then
  scatter_df_x <- scatter_df_x[order(rownames(scatter_df_x)),]
  scatter_df_x <- scatter_df_x[which(rownames(scatter_df_x) %in% combined_x$gene),]
  combined_x <- combined_x[which(combined_x$gene %in% rownames(scatter_df_x)),] 
  #Bind autosomal significance columns to scatter_dfs
  scatter_df <- cbind(scatter_df, combined[which(combined$gene %in% rownames(scatter_df)),c("auto_sig", "p_score", "auto_beta")])
  scatter_df_x <- cbind(scatter_df_x, combined_x[which(combined_x$gene %in% rownames(scatter_df_x)),c("auto_sig", "p_score", "auto_beta")])
  
  
  #Make volcano plots
  #Autosomal first
  tiff(paste(out_dir, "volcano_plots/", gender, "_auto_volcano.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  suppressWarnings(print(ggplot(combined[which(combined$Chr_type == "Autosomal"),], aes(x = get(paste("auto_", gender, "_FC", sep="")), y=neg_log_auto)) + #volcanoplot with log2Foldchange versus pvalue
                           geom_point(aes(col=auto_sig)) + #add points colored by significance
                           scale_color_manual(values=c("red", "black")) + 
                           ylim(0,18) + 
                           xlim(-4,4) + 
                           theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                           xlab("logFC") + ylab("-log adj. P-Value") + 
                           ggtitle(paste(capitalize(gender), " DE Genes in Response to LPS", sep = "")))) #e.g. 'Volcanoplot DESeq2')
  dev.off()
  
  #X second from autosomal
  tiff(paste(out_dir, "volcano_plots/", gender, "_x_volcano.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  suppressWarnings(print(ggplot(combined[which(combined$Chr_type == "X-linked"),], aes(x = get(paste("auto_", gender, "_FC", sep="")), y=neg_log_auto)) + #volcanoplot with log2Foldchange versus pvalue
                           geom_point(aes(col=auto_sig)) + #add points colored by significance
                           scale_color_manual(values=c("red", "black")) + 
                           ylim(0,18) + 
                           xlim(-4,4) + 
                           theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                           xlab("logFC") + ylab("-log adj. P-Value") + 
                           ggtitle(paste(capitalize(gender), " DE Genes in Response to LPS", sep = "")))) #e.g. 'Volcanoplot DESeq2')
  dev.off()
  
  #X third from X-GRM
  tiff(paste(out_dir, "volcano_plots/", gender, "_x_volcano_x_grm.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  suppressWarnings(print(ggplot(combined_x[which(combined_x$Chr_type == "X-linked"),], aes(x = get(paste("auto_", gender, "_FC", sep="")), y=neg_log_auto)) + #volcanoplot with log2Foldchange versus pvalue
                           geom_point(aes(col=auto_sig)) + #add points colored by significance
                           scale_color_manual(values=c("red", "black")) + 
                           ylim(0,18) + 
                           xlim(-4,4) + 
                           theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                           xlab("logFC") + ylab("-log adj. P-Value") + 
                           ggtitle(paste(capitalize(gender), " DE Genes in Response to LPS", sep = "")))) #e.g. 'Volcanoplot DESeq2')
  dev.off()
  
  #Match data frame rows and sort
  pvals<-pvals[which(pvals$gene %in% rownames(scatter_df)),]
  pvals_x<-pvals_x[which(pvals_x$gene %in% rownames(scatter_df)),]
  scatter_df<-scatter_df[which(rownames(scatter_df) %in% pvals$gene),]
  pvals<-pvals[order(pvals$gene),]
  pvals_x<-pvals_x[order(pvals_x$gene),]
  scatter_df<-scatter_df[order(rownames(scatter_df)),]
  #Create data frame of x-genes
  x_specific <- cbind.data.frame(gene=pvals$gene, auto_p=pvals$p_score, x_p=pvals_x$p_score, sig_auto=pvals$sig, sig_x=pvals_x$sig, Chr_type=scatter_df$Chr_type)
  x_specific <- x_specific[which(x_specific$Chr_type == "X-linked"),] 
  x_specific <- cbind(x_specific, neg_log_x=-log(x_specific[,"x_p"], base = 10))
  x_specific <- cbind(x_specific, neg_log_auto=-log(x_specific[,"auto_p"], base = 10))
  x_specific <- cbind(x_specific, significance=ifelse(x_specific$sig_x == "not significant" & x_specific$sig_auto == "not significant", "Neither", 
                                                          ifelse(x_specific$sig_x != "not significant" & x_specific$sig_auto == "not significant", "X-GRM Only", 
                                                                 ifelse(x_specific$sig_x == "not significant" & x_specific$sig_auto != "not significant", "Auto-GRM Only", "Both GRMs"))))
  
  
  #Make plots
  tiff(paste(out_dir, "comparison_scatter_plots/", gender, "_de_comp.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  print(suppressWarnings(ggplot(x_specific, aes(x = auto_p, y=x_p)) + 
                           geom_point(aes(col=significance)) + #add points colored by significance
                           scale_color_manual(values=c("#E95C20FF", "#4F2C1DFF", "lightgray", "#00ff00")) + 
                           ylim(0,1) + 
                           xlim(0,1) + 
                           theme_classic() + 
                           geom_abline(intercept = 0, slope = 1, color = "lightgray", linetype="dashed", size=1.2) + 
                           #geom_segment(aes(x = -1, y = -1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = -1, y = -1, xend = 1, yend = -1), linetype="dashed", color="lightgray") + 
                           #geom_segment(aes(x = 1, y = 1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = 1, y = 1, xend = 1, yend = -1), linetype="dashed", color="lightgray") +
                           theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                           #geom_text_repel(aes(label=gene, color = factor(pvals$sig)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                           #geom_hline(yintercept = 0, linetype="dashed", color="lightgray") + geom_vline(xintercept = 0, linetype="dashed", color="lightgray") + 
                           #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                           xlab("FDR Adj. P-Value (Auto-GRM)") + ylab("FDR Adj. P-Value (X-GRM)")))
  dev.off()
  tiff(paste(out_dir, "comparison_scatter_plots/", gender, "_de_comp_zoom.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  suppressWarnings(print(ggplot(x_specific, aes(x = auto_p, y=x_p)) + 
                           geom_point(aes(col=significance)) + #add points colored by significance
                           scale_color_manual(values=c("#E95C20FF", "#4F2C1DFF", "lightgray", "#00ff00")) + 
                           ylim(0,0.15) + 
                           xlim(0,0.15) + 
                           theme_classic() + 
                           geom_abline(intercept = 0, slope = 1, color = "lightgray", linetype="dashed", size=1.2) + 
                           #geom_segment(aes(x = -1, y = -1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = -1, y = -1, xend = 1, yend = -1), linetype="dashed", color="lightgray") + 
                           #geom_segment(aes(x = 1, y = 1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = 1, y = 1, xend = 1, yend = -1), linetype="dashed", color="lightgray") +
                           theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                           #geom_text_repel(aes(label=gene, color = factor(pvals$sig)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                           #geom_hline(yintercept = 0, linetype="dashed", color="lightgray") + geom_vline(xintercept = 0, linetype="dashed", color="lightgray") + 
                           #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                           xlab("FDR Adj. P-Value (Auto-GRM)") + ylab("FDR Adj. P-Value (X-GRM)")))
  dev.off()
  tiff(paste(out_dir, "comparison_scatter_plots/", gender, "_de_comp_log.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  print(suppressWarnings(ggplot(x_specific, aes(x = neg_log_auto, y=neg_log_x)) + 
                           geom_point(aes(col=significance)) + #add points colored by significance
                           scale_color_manual(values=c("#E95C20FF", "#4F2C1DFF", "lightgray", "#00ff00")) + 
                           ylim(0,max(x_specific$neg_log_x)) + 
                           xlim(0,max(x_specific$neg_log_auto)) + 
                           theme_classic() + 
                           geom_abline(intercept = 0, slope = 1, color = "lightgray", linetype="dashed", size=1.2) + 
                           #geom_segment(aes(x = -1, y = -1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = -1, y = -1, xend = 1, yend = -1), linetype="dashed", color="lightgray") + 
                           #geom_segment(aes(x = 1, y = 1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = 1, y = 1, xend = 1, yend = -1), linetype="dashed", color="lightgray") +
                           theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                           #geom_text_repel(aes(label=gene, color = factor(pvals$sig)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                           #geom_hline(yintercept = 0, linetype="dashed", color="lightgray") + geom_vline(xintercept = 0, linetype="dashed", color="lightgray") + 
                           #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                           xlab("-log(FDR Adj. P-Value (Auto-GRM))") + ylab("-log(FDR Adj. P-Value (X-GRM)")))
  dev.off()

  #Make pdf for Michelle
  pdf(file = paste(out_dir, "comparison_scatter_plots/", "Supplemental_Figure_4", ifelse(gender == "females", "B", "A"), ".pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
  print(suppressWarnings(ggplot(x_specific, aes(x = neg_log_auto, y=neg_log_x)) + 
                           geom_point(aes(col=significance)) + #add points colored by significance
                           scale_color_manual(values=c("#E95C20FF", "#4F2C1DFF", "lightgray", "#00ff00")) + 
                           ylim(0,max(x_specific$neg_log_x)) + 
                           xlim(0,max(x_specific$neg_log_auto)) + 
                           theme_classic() + 
                           geom_abline(intercept = 0, slope = 1, color = "lightgray", linetype="dashed", size=1.2) + 
                           #geom_segment(aes(x = -1, y = -1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = -1, y = -1, xend = 1, yend = -1), linetype="dashed", color="lightgray") + 
                           #geom_segment(aes(x = 1, y = 1, xend = -1, yend = 1), linetype="dashed", color="lightgray") + geom_segment(aes(x = 1, y = 1, xend = 1, yend = -1), linetype="dashed", color="lightgray") +
                           theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                           #geom_text_repel(aes(label=gene, color = factor(pvals$sig)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                           #geom_hline(yintercept = 0, linetype="dashed", color="lightgray") + geom_vline(xintercept = 0, linetype="dashed", color="lightgray") + 
                           annotate("text", x=min(x_specific$neg_log_auto) + 0.5, y=max(x_specific$neg_log_x), label = paste("~R^2==~", formatC(summary(lm(x_specific$neg_log_x ~ x_specific$neg_log_auto))$r.squared, digits=3)), parse = TRUE) +
                           xlab("-log(FDR Adj. P-Value (Auto-GRM))") + ylab("-log(FDR Adj. P-Value (X-GRM)")))
  dev.off()
  
  #Create export of gene lists
  export_list <- cbind(as.character(pvals$gene), as.character(pvals$sig), as.numeric(pvals$p_score))
  export_list2 <- cbind(as.character(pvals_x$gene), as.character(pvals_x$sig), as.numeric(pvals_x$p_score))
  table(export_list[,2])
  colnames(export_list) <- c("gene", "sig", "p_score")
  colnames(export_list2) <- c("gene", "sig", "p_score")
  #sort export list
  export_list <- export_list[order(export_list[,2]),]
  export_list2 <- export_list2[order(export_list2[,2]),]
  #Write to file
  write.table(export_list, file = paste("C:/Users/Janel/Documents/sex-specific/de_analysis/", gender, "_sig_global_correction.tsv", sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(export_list2, file = paste("C:/Users/Janel/Documents/sex-specific/de_analysis/", gender, "_sig_global_correction_x_grm.tsv", sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}



# 4) Comparison of fold-changes ====

#Remove junk data
suppressWarnings(remove(auto_raw, combined, export_list, foldchange_df, individuals_bound, pvals, foldchange_df_x, foldchange_df_auto, auto_raw_auto))

#Rename columns of scatter_df
colnames(scatter_df) <- c("females_FC", "males_FC", "Chr_type", "females_sig", "females_p_score", "females_auto_beta", "males_sig", "males_p_score", "males_auto_beta")

#Make joint significance column
scatter_df_add <- cbind(scatter_df, significance=ifelse(scatter_df$females_sig == "Not Sig" & scatter_df$males_sig == "Not Sig", "Neither", 
                                                        ifelse(scatter_df$females_sig != "Not Sig" & scatter_df$males_sig == "Not Sig", "Females Only", 
                                                               ifelse(scatter_df$females_sig == "Not Sig" & scatter_df$males_sig != "Not Sig", "Males Only", "Both Sexes"))))


#Create color-coded FC Scatter Data and add labels to data frame
input_add<-cbind(scatter_df_add, Expression=ifelse(abs(scatter_df_add[,"males_FC"]) > abs(scatter_df_add[,"females_FC"]), "Males", "Females"))

#Write FC and p-values for Emma
write_file <- cbind.data.frame(gene=rownames(input_add), Chr_type= input_add$Chr_type, significance=input_add$significance, 
                               females_FC=input_add$females_FC, males_FC=input_add$males_FC,
                               females_p_score=input_add$females_p_score, males_p_score=input_add$males_p_score)
write.table(write_file, file = "C:/Users/Janel/Documents/sex-specific/de_analysis/foldchanges.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#Plot FC data
tiff(paste(out_dir, "fold_change_scatter_plots/", "Fold_Change_Comp_FC.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_add[which(scatter_df_add$Chr_type == "X-linked"),], aes(x = females_FC, y = males_FC)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(-8,12) + 
                   xlim(-8,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (fold change)") + ylab("Male response to LPS \n (fold change)"))
dev.off()
tiff(paste(out_dir, "fold_change_scatter_plots/", "Fold_Change_Comp_autosome_FC.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_add[which(scatter_df_add$Chr_type == "Autosomal"),], aes(x = females_FC, y = males_FC)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(-8,12) + 
                   xlim(-8,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (fold change)") + ylab("Male response to LPS \n (fold change)"))
dev.off()
#Make PDF versions for Michelle
pdf(paste(out_dir, "fold_change_scatter_plots/", "Figure_3B.pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
suppressWarnings(ggplot(scatter_df_add[which(scatter_df_add$Chr_type == "X-linked"),], aes(x = females_FC, y = males_FC)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(-8,12) + 
                   xlim(-8,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (fold change)") + ylab("Male response to LPS \n (fold change)"))
dev.off()
pdf(paste(out_dir, "fold_change_scatter_plots/", "Figure_3A.pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
suppressWarnings(ggplot(scatter_df_add[which(scatter_df_add$Chr_type == "Autosomal"),], aes(x = females_FC, y = males_FC)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(-8,12) + 
                   xlim(-8,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (fold change)") + ylab("Male response to LPS \n (fold change)"))
dev.off()


# 5) Make Plots of Absolute Values of Fold Change ====

#Bind abs value of fold change columns onto scatter_df_add
scatter_df_add_two <- cbind(scatter_df_add, females_FC_abs=abs(scatter_df_add$females_FC), males_FC_abs=abs(scatter_df_add$males_FC))

#Calculate regression lines for plot
complete_lm<- lm(scatter_df_add_two[,"males_FC_abs"] ~ scatter_df_add_two[,"females_FC_abs"] - 1)
auto_lm <- lm(scatter_df_add_two[which(scatter_df_add_two$Chr_type == "Autosomal"),"males_FC_abs"] ~ scatter_df_add_two[which(scatter_df_add_two$Chr_type == "Autosomal"),"females_FC_abs"] - 1)
x_lm <- lm(scatter_df_add_two[which(scatter_df_add_two$Chr_type == "X-linked"),"males_FC_abs"] ~ scatter_df_add_two[which(scatter_df_add_two$Chr_type == "X-linked"),"females_FC_abs"] - 1)

#Plot abs FC data
tiff(paste(out_dir, "abs_fold_change_scatter_plots/", "Fold_Change_Comp_abs_all_genes_FC.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_add_two, aes(x = females_FC_abs, y = males_FC_abs)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(0,12) + 
                   xlim(0,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   geom_abline(slope = complete_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("black", 0.4)) + #all genes regression line
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (absolute value of fold change)") + ylab("Male response to LPS \n (absolute value of fold change)"))
dev.off()
tiff(paste(out_dir, "abs_fold_change_scatter_plots/", "Fold_Change_Comp_abs_x-linked_FC.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_add_two[which(scatter_df_add$Chr_type == "X-linked"),], aes(x = females_FC_abs, y = males_FC_abs)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(0,12) + 
                   xlim(0,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   geom_abline(slope = x_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("black", 0.4)) + #x regression line
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (absolute value of fold change)") + ylab("Male response to LPS \n (absolute value of fold change)"))
dev.off()
tiff(paste(out_dir, "abs_fold_change_scatter_plots/", "Fold_Change_Comp_abs_autosome_FC.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_add_two[which(scatter_df_add$Chr_type == "Autosomal"),], aes(x = females_FC_abs, y = males_FC_abs)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(0,12) + 
                   xlim(0,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   geom_abline(slope = auto_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("black", 0.4)) + #autosomal regression line
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (absolute value of fold change)") + ylab("Male response to LPS \n (absolute value of fold change)"))
dev.off()

# 5.1) Make graphs with X-linked influential points removed and marked ====

#Label points in df with abs fold changes greater than 3.0 in both dimensions (males and females)
scatter_df_add_three <- cbind.data.frame(scatter_df_add_two, gene_labels=ifelse(scatter_df_add_two$females_FC_abs >= 3 & scatter_df_add_two$males_FC_abs >= 3, rownames(scatter_df_add_two), ""))

#Plot graph
tiff(paste(out_dir, "abs_fold_change_scatter_plots/", "Fold_Change_Comp_abs_x-linked_FC_labeled.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_add_three[which(scatter_df_add_three$Chr_type == "X-linked"),], aes(x = females_FC_abs, y = males_FC_abs)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(0,12) + 
                   xlim(0,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   geom_text_repel(aes(label=gene_labels, color = significance), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   geom_abline(slope = x_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("black", 0.4)) + #x regression line
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (absolute value of fold change)") + ylab("Male response to LPS \n (absolute value of fold change)"))
dev.off()

#Remove labeled points from scatter_df_add_three and regraph
scatter_df_add_three <- scatter_df_add_three[which(scatter_df_add_three$gene_labels == ""),]
#Recalculate regression line
new_x_lm <- lm(scatter_df_add_three[which(scatter_df_add_three$Chr_type == "X-linked"),"males_FC_abs"] ~ scatter_df_add_three[which(scatter_df_add_three$Chr_type == "X-linked"),"females_FC_abs"] - 1)

#Plot graph
#Plot graph
tiff(paste(out_dir, "abs_fold_change_scatter_plots/", "Fold_Change_Comp_abs_x-linked_FC_no_outliers.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_add_three[which(scatter_df_add_three$Chr_type == "X-linked"),], aes(x = females_FC_abs, y = males_FC_abs)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(0,12) + 
                   xlim(0,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=gene_labels, color = significance), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   geom_abline(slope = new_x_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("black", 0.4)) + #x regression line
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (absolute value of fold change)") + ylab("Male response to LPS \n (absolute value of fold change)"))
dev.off()

#Delete scatter_df_add_three
remove(scatter_df_add_three)

# 5.2) Look at escape genes====

#Read in escape genes
blood_escape <- read.table("C:/Users/Janel/Documents/sex-specific/de_analysis/WHLBLOOD_escape_genes_121619.txt", header = FALSE, stringsAsFactors = FALSE)
#Create and filter scatter_df_escape
scatter_df_escape <- cbind.data.frame(scatter_df_add_two, escape_status=ifelse(rownames(scatter_df_add_two) %in% blood_escape[,1], TRUE, FALSE))
scatter_df_escape <- scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked"),]

#Calculate regression lines for plot
x_escape_lm <- lm(scatter_df_escape[which(scatter_df_escape$escape_status == TRUE),"males_FC_abs"] ~ scatter_df_escape[which(scatter_df_escape$escape_status == TRUE),"females_FC_abs"] - 1)
x_inactivated_lm <- lm(scatter_df_escape[which(scatter_df_escape$escape_status == FALSE),"males_FC_abs"] ~ scatter_df_escape[which(scatter_df_escape$escape_status == FALSE),"females_FC_abs"] - 1)

#Make graphs
tiff(paste(out_dir, "abs_fold_change_scatter_plots/", "Fold_Change_abs_x-linked_escape_FC.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_escape[which(scatter_df_escape$escape_status == TRUE),], aes(x = females_FC_abs, y = males_FC_abs)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(0,12) + 
                   xlim(0,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   geom_abline(slope = x_escape_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("black", 0.4)) + #x regression line
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (absolute value of fold change)") + ylab("Male response to LPS \n (absolute value of fold change)"))
dev.off()
tiff(paste(out_dir, "abs_fold_change_scatter_plots/", "Fold_Change_abs_x-linked_inactivated_FC.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
suppressWarnings(ggplot(scatter_df_escape[which(scatter_df_escape$escape_status == FALSE),], aes(x = females_FC_abs, y = males_FC_abs)) + 
                   geom_point(aes(col=significance)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("purple", "red", "blue", "lightgray")) +
                   ylim(0,12) + 
                   xlim(0,12) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1) + 
                   geom_abline(slope = x_inactivated_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("black", 0.4)) + #x regression line
                   #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                   xlab("Female response to LPS \n (absolute value of fold change)") + ylab("Male response to LPS \n (absolute value of fold change)"))
dev.off()

# 6) Wilcoxon tests on Fold changes ====

#Remove some junk
suppressWarnings(remove(scatter_df_add, scatter_df_escape, write_file, x_escape_lm, x_inactivated_lm, x_lm, x_specific, new_x_lm, complete_lm, auto_lm))
suppressWarnings(remove(pvals_x, input_add, export_list2, combined_x, auto_raw_x, gender))

#Remake scatter_df_escape but with all genes
scatter_df_escape <- cbind.data.frame(scatter_df_add_two, escape_status=ifelse(rownames(scatter_df_add_two) %in% blood_escape[,1], TRUE, FALSE))
remove(scatter_df_add_two)

#Run tests
all_test <- wilcox.test(scatter_df_escape[,"females_FC_abs"], scatter_df_escape[,"males_FC_abs"], paired = TRUE, alternative = "two.sided")
autosomal_test <- wilcox.test(scatter_df_escape[which(scatter_df_escape$Chr_type == "Autosomal"),"females_FC_abs"], 
                              scatter_df_escape[which(scatter_df_escape$Chr_type == "Autosomal"),"males_FC_abs"], paired = TRUE, 
                              alternative = "two.sided")
x_test <- wilcox.test(scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked"),"females_FC_abs"], 
                              scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked"),"males_FC_abs"], paired = TRUE, 
                              alternative = "two.sided")
escape_test <- wilcox.test(scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked" & scatter_df_escape$escape_status == TRUE),"females_FC_abs"], 
                           scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked" & scatter_df_escape$escape_status == TRUE),"males_FC_abs"], 
                           paired = TRUE, alternative = "two.sided")
inactivated_test <- wilcox.test(scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked" & scatter_df_escape$escape_status == FALSE),"females_FC_abs"], 
                           scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked" & scatter_df_escape$escape_status == FALSE),"males_FC_abs"], 
                           paired = TRUE, alternative = "two.sided")
autosomal_de_both_test <- wilcox.test(scatter_df_escape[which(scatter_df_escape$Chr_type == "Autosomal" & scatter_df_escape$significance == "Both Sexes"),"females_FC_abs"], 
                                      scatter_df_escape[which(scatter_df_escape$Chr_type == "Autosomal" & scatter_df_escape$significance == "Both Sexes"),"males_FC_abs"], 
                                      paired = TRUE, alternative = "two.sided")

# 7) Calculate "True" Sex-Specificity ====  

#Calculate True Sex-specificity and append columns onto scatter_df_escape
scatter_df_escape <- cbind.data.frame(scatter_df_escape, 
                                      true_sig=ifelse(scatter_df_escape$males_sig != "Not Sig" & scatter_df_escape$females_p_score > 0.2, "Male-Specific",
                                                      ifelse(scatter_df_escape$females_sig != "Not Sig" & scatter_df_escape$males_p_score > 0.2, "Female-Specific",
                                                             "Not Specific")))
#Count males and females true sig results
x_specific_results <- sum(nrow(scatter_df_escape[which(scatter_df_escape$true_sig == "Male-Specific" & scatter_df_escape$Chr_type == "X-linked"),]) + nrow(scatter_df_escape[which(scatter_df_escape$true_sig == "Female-Specific" & scatter_df_escape$Chr_type == "X-linked"),]))
x_non_specific_results <- sum(nrow(scatter_df_escape[which(scatter_df_escape$true_sig == "Not Specific" & scatter_df_escape$Chr_type == "X-linked"),]))
auto_specific_results <- sum(nrow(scatter_df_escape[which(scatter_df_escape$true_sig == "Male-Specific" & scatter_df_escape$Chr_type == "Autosomal"),]) + nrow(scatter_df_escape[which(scatter_df_escape$true_sig == "Female-Specific" & scatter_df_escape$Chr_type == "Autosomal"),]))
auto_non_specific_results <- sum(nrow(scatter_df_escape[which(scatter_df_escape$true_sig == "Not Specific" & scatter_df_escape$Chr_type == "Autosomal"),]))

#Run test and return results
sex_specific_table <- matrix(c(x_specific_results, x_non_specific_results, auto_specific_results, auto_non_specific_results), byrow = FALSE, nrow = 2, ncol = 2, dimnames = list(c("males", "females"), c("X-linked", "Autosomal")))
#Run test
sex_specific_test <- fisher.test(sex_specific_table) 
#Clean-Up
remove(x_specific_results, x_non_specific_results, auto_specific_results, auto_non_specific_results)

#Run escape test for sex_specific results
sex_specific_escape_table <- matrix(c(nrow(scatter_df_escape[which(scatter_df_escape$escape_status == TRUE & scatter_df_escape$true_sig != "Not Specific"),]), 
                                      nrow(scatter_df_escape[which(scatter_df_escape$escape_status == TRUE & scatter_df_escape$true_sig == "Not Specific"),]), 
                                      nrow(scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked" & scatter_df_escape$true_sig != "Not Specific"),]),
                                      nrow(scatter_df_escape[which(scatter_df_escape$Chr_type == "X-linked" & scatter_df_escape$true_sig == "Not Specific"),])), byrow = FALSE, nrow = 2, ncol = 2, dimnames = list(c("sex-specific", "non-specific"), c("Escape", "All X-linked")))
sex_specific_escape_test <- fisher.test(sex_specific_escape_table)

#Write sex-specific lists to file for supplemental table 4
write.table(rownames(scatter_df_escape[which(scatter_df_escape$true_sig == "Female-Specific" & scatter_df_escape$Chr_type == "X-linked"),]), file = paste(inp_dir, "Female-Specific_Genes.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(scatter_df_escape[which(scatter_df_escape$true_sig == "Male-Specific" & scatter_df_escape$Chr_type == "X-linked"),]), file = paste(inp_dir, "Male-Specific_Genes.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(scatter_df_escape[which(scatter_df_escape$significance != "Both Sexes" & 
                                               scatter_df_escape$Chr_type == "X-linked" &
                                               scatter_df_escape$females_p_score <= 0.2 & 
                                               scatter_df_escape$males_p_score <= 0.2),]), 
            file = paste(inp_dir, "twenty_pct_other.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 8) Run Sex enrichment tests ====

#Create quick function for testing all sig results, so that test can be duplicated in all sets
run_proportion_test <- function(data_frame, Chr_type){

  #Count males and females sig results
  males_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == Chr_type), "females_sig"], data_frame[which(data_frame$Chr_type == Chr_type), "males_sig"])[,1])
  females_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == Chr_type), "females_sig"], data_frame[which(data_frame$Chr_type == Chr_type), "males_sig"])[1,])
  sig_results <-  c(males_sig_results, females_sig_results)
  
  #Create proportions
  proportion_denominator <- c(nrow(data_frame[which(data_frame$Chr_type == Chr_type),]), nrow(data_frame[which(data_frame$Chr_type == Chr_type),])) 

  #Run test and return results
  return(prop.test(sig_results, proportion_denominator, alternative = "two.sided", correct = TRUE))
}

#Create quick function for testing only sex-specific results, so that test can be duplicated in all sets
run_specific_proportion_test <- function(data_frame, Chr_type){
  
  #Count males and females sig results
  males_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == Chr_type), "females_sig"], data_frame[which(data_frame$Chr_type == Chr_type), "males_sig"])[2,1])
  females_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == Chr_type), "females_sig"], data_frame[which(data_frame$Chr_type == Chr_type), "males_sig"])[1,2])
  sig_results <-  c(males_sig_results, females_sig_results)
  
  #Create proportions
  proportion_denominator <- c(nrow(data_frame[which(data_frame$Chr_type == Chr_type),]), nrow(data_frame[which(data_frame$Chr_type == Chr_type),])) 
  
  #Run test and return results
  return(prop.test(sig_results, proportion_denominator, alternative = "two.sided", correct = TRUE))
}

#Create quick function for testing only sex-specific results, so that test can be duplicated in all sets
run_enrichment_test <- function(data_frame){
  
  #Count males and females sig results
  males_x_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == "X-linked"), "females_sig"], data_frame[which(data_frame$Chr_type == "X-linked"), "males_sig"])[,1])
  females_x_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == "X-linked"), "females_sig"], data_frame[which(data_frame$Chr_type == "X-linked"), "males_sig"])[1,])
  males_auto_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == "Autosomal"), "females_sig"], data_frame[which(data_frame$Chr_type == "Autosomal"), "males_sig"])[,1])
  females_auto_sig_results <- sum(table(data_frame[which(data_frame$Chr_type == "Autosomal"), "females_sig"], data_frame[which(data_frame$Chr_type == "Autosomal"), "males_sig"])[1,])
  
  #Run test and return results
  return(fisher.test(matrix(c(males_auto_sig_results, females_auto_sig_results, males_x_sig_results, females_x_sig_results), byrow = FALSE, nrow = 2, ncol = 2, dimnames = list(c("males", "females"), c("Autosomal", "X-linked")))))
}

#Run tests
autosomal_sex_test <- run_proportion_test(scatter_df_escape, "Autosomal")
x_sex_test <- run_proportion_test(scatter_df_escape, "X-linked")
autosomal_sex_specific_test <- run_specific_proportion_test(scatter_df_escape, "Autosomal")
x_sex_specific_test <- run_specific_proportion_test(scatter_df_escape, "X-linked")
sex_enrichment <- run_enrichment_test(scatter_df_escape)


# 9) Check directions ====

#Remove junk
suppressWarnings(remove(input_add, filter_dir, gender))

#make tables to check de categories
with(scatter_df, table(females_sig, males_sig))

#bind on significance column
scatter_df <- cbind(scatter_df, significance=ifelse(scatter_df$females_sig == "Not Sig" & scatter_df$males_sig == "Not Sig", "Neither", 
                                                    ifelse(scatter_df$females_sig != "Not Sig" & scatter_df$males_sig == "Not Sig", "Females Only", 
                                                           ifelse(scatter_df$females_sig == "Not Sig" & scatter_df$males_sig != "Not Sig", "Males Only", "Both Sexes"))))
#Reverse direction of male betas to account for Gemma error
scatter_df$males_auto_beta <- -scatter_df$males_auto_beta
scatter_df_escape$males_auto_beta <- -scatter_df_escape$males_auto_beta

#identify genes with opposite directions
opposite_df <- scatter_df[which(abs(scatter_df$females_auto_beta+scatter_df$males_auto_beta)<abs(scatter_df$females_auto_beta)+abs(scatter_df$males_auto_beta)),]
#make table for opposite df
table(opposite_df$significance)

#Check for DE results with opposite directions
opposite_DE_df <- opposite_df[which(opposite_df$Chr_type == "X-linked" & opposite_df$significance != "Neither"),]
#Append escape info
opposite_DE_df <- cbind.data.frame(opposite_DE_df, escape_status=ifelse(rownames(opposite_DE_df) %in% blood_escape[,1], TRUE, FALSE))

#Calculate abs fold change stat
pct_male_abs_change <- nrow(scatter_df[which(abs(scatter_df$males_FC)>abs(scatter_df$females_FC) & scatter_df$significance == "Both Sexes"),])/nrow(scatter_df[which(scatter_df$significance == "Both Sexes"),])
pct_male_abs_change_auto <- nrow(scatter_df[which(abs(scatter_df$males_FC)>abs(scatter_df$females_FC) & scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "Autosomal"),])/nrow(scatter_df[which(scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "Autosomal"),])
pct_male_abs_change_x <- nrow(scatter_df[which(abs(scatter_df$males_FC)>abs(scatter_df$females_FC) & scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "X-linked"),])/nrow(scatter_df[which(scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "X-linked"),])

#Run fisher's exact test for x-linked and autosomal male FC > female FC
#make table for test
pct_male_table <- matrix(c(nrow(scatter_df[which(abs(scatter_df$males_FC)>abs(scatter_df$females_FC) & scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "Autosomal"),])
                        ,nrow(scatter_df[which(abs(scatter_df$males_FC)<abs(scatter_df$females_FC) & scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "Autosomal"),]),
                        nrow(scatter_df[which(abs(scatter_df$males_FC)>abs(scatter_df$females_FC) & scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "X-linked"),]),
                        nrow(scatter_df[which(abs(scatter_df$males_FC)<abs(scatter_df$females_FC) & scatter_df$significance == "Both Sexes" & scatter_df$Chr_type == "X-linked"),])), 
                        nrow = 2, ncol = 2,byrow = FALSE)
colnames(pct_male_table) <- c("Autosomal", "X-linked")
rownames(pct_male_table) <- c("males>females", "females>males")
#Execute tests
pct_male_test <- fisher.test(pct_male_table)

# 10) Read in x-chr results for Fisher's exact tests ====

#Clear junk
suppressWarnings(remove(pvals, genders))

#Bind scatter table together
scatter_table <- table(scatter_df_escape$true_sig, scatter_df_escape$Chr_type)
colnames(scatter_table) <- c("X-Chr", "Autosomes")

#Create binary_table
binary_table <- rbind(scatter_table[3,], scatter_table[1,] + scatter_table[2,])
rownames(binary_table) <- c("Non-Specific", "Gender-Specific")

#Perform Fisher's exact tests
binary_test <- fisher.test(binary_table)
scatter_test <- fisher.test(scatter_table)

# 11) Print out gene lists for X-chromosome based on gender signficance ====

write.table(rownames(scatter_df_escape[which(scatter_df_escape$significance == "Both Sexes" & scatter_df_escape$Chr_type == "X-linked"),]), 
            file = paste(inp_dir, "X_Linked_Both_Sexes.txt", sep = ""), quote = FALSE, col.names = FALSE, row.names = FALSE)


# 12) Examine HLA genes in particular ====

#Read in expr file and get list of expressed HLA genes
expressed_genes <- read.table("C:/Users/Janel/Documents/sex-specific/processed_data/BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", header = TRUE, stringsAsFactors = FALSE)
expressed_genes <- rownames(expressed_genes)
#Filter out junk genes
expressed_genes <- expressed_genes[which(expressed_genes %in% junk_genes == FALSE)]
#Filter to HLA_genes
expressed_HLA <- expressed_genes[which(regexpr("HLA-", expressed_genes) > 0)]

#Make scatter_df_HLA
scatter_df_HLA <- scatter_df_escape[expressed_HLA,]
#There are no sex-specific DE HLA genes.  There also don't appear to be any specific patterns of what is and is not DE.  
#DE genes bridge class I and II as well as classical and non-classical

#Check for enrichment of DE genes among HLA genes relative to all autosomal genes
HLA_test_table <- cbind.data.frame(HLA = c(length(expressed_HLA), nrow(scatter_df_HLA[which(scatter_df_HLA$significance != "Neither"),])), 
                                   ALL = c(nrow(scatter_df_escape[which(scatter_df_escape$Chr_type =="Autosomal"),]), 
                                           nrow(scatter_df_escape[which(scatter_df_escape$Chr_type =="Autosomal" & 
                                                                          scatter_df_escape$significance !="Neither"),])))
HLA_test <- fisher.test(HLA_test_table)
#No significant enrichment or depletion of DE results among HLA genes relative to all autosomal genes
#Positive Betas are upregulated
                                   
#Clear some junk
suppressWarnings(remove(opposite_DE_df, opposite_df, scatter_df, pheno_HLA, scatter_df_x, HLA_genes, inactivated_test, scatter_table, scatter_test))

#Read in sex-by-treatment data to get info on 1055 sbt genes     
sex_by_treatment <- read.table(file = "C:/Users/Janel/Documents/sex-specific/overall_analyses/sbt.all.gemma.assoc.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE, row.names = 1)
colnames(sex_by_treatment) <- c("chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "p_score")
#Remove junk genes from sex_by_treatment
sex_by_treatment <- sex_by_treatment[which(rownames(sex_by_treatment) %in% junk_genes == FALSE),]
#Append p_adj_score column
sex_by_treatment <- cbind.data.frame(sex_by_treatment, p_adj_score=p.adjust(sex_by_treatment$p_score, "fdr"))
#Append p_adj_sbt to scatter_df_HLA
scatter_df_HLA <- cbind.data.frame(scatter_df_HLA, p_adj_sbt=sex_by_treatment[rownames(scatter_df_HLA),"p_adj_score"])
#Clear junk
remove(sex_by_treatment)
#Append two columns one that Checks for opposite directionality of fold changes and one that compares abs fold change
scatter_df_HLA <- cbind.data.frame(scatter_df_HLA,
                                   opposite=ifelse(abs(scatter_df_HLA$females_FC) + abs(scatter_df_HLA$males_FC) > 
                                                     abs(scatter_df_HLA$females_FC + scatter_df_HLA$males_FC), TRUE, FALSE), 
                                   abs_FC_comp=ifelse(scatter_df_HLA$females_FC_abs < scatter_df_HLA$males_FC, 
                                                      "males greater", "females greater"))

#Make an export file of HLA genes
write.table(scatter_df_HLA[which(scatter_df_HLA$p_adj_sbt < 0.05),c("females_FC", "males_FC", "opposite", "abs_FC_comp")], 
            file = paste(inp_dir, "HLA_sbt_table.txt", sep = ""), col.names = TRUE, row.names = TRUE, quote = FALSE)




