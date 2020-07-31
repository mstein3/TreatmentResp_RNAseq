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


# 1) Read in and Analyze Overall Results ====

#Set directories
inp_dir <- "C:/Users/Janel/Documents/sex-specific/overall_analyses/"
filter_dir <- "C:/Users/Janel/Documents/sex-specific/filters/"
#Set junk genes
junk_genes <- c("HLA-DPA1", "HLA-DRB5", "HLA-G")

#Read in results
overall_raw <- read.table(paste(inp_dir, "overall.all.gemma.assoc.txt", sep = ""), sep = "\t", header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(overall_raw) <- c("chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score")
#Remove junk_genes
overall_raw <- overall_raw[which(rownames(overall_raw) %in% junk_genes == FALSE),]

#Bind on adj_p_score column
overall_raw <- cbind.data.frame(overall_raw, adj_p_score=p.adjust(overall_raw$p_score, "fdr"))

#Calculate significant results
num_overall <- nrow(overall_raw[which(overall_raw$adj_p_score < 0.05),])

# 2) Read in and Sex by Treatment Results ====

#Read in results
sbt_raw <- read.table(paste(inp_dir, "sbt.all.gemma.assoc.txt", sep = ""), sep = "\t", header = FALSE, row.names = 1, stringsAsFactors = FALSE)
colnames(sbt_raw) <- c("chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "logl_H1", "p_score")
#Remove junk_genes
sbt_raw <- sbt_raw[which(rownames(sbt_raw) %in% junk_genes == FALSE),]

#Bind on adj_p_score column
sbt_raw <- cbind.data.frame(sbt_raw, adj_p_score=p.adjust(sbt_raw$p_score, "fdr"))

#Calculate significant results
num_sbt <- nrow(sbt_raw[which(sbt_raw$adj_p_score < 0.05),])

#Read in Michelle's results and compare to mine
michelle_raw <- read.table(paste(inp_dir, "CombinedResults_TreatbySexMainEffect_wFDR.txt", sep = ""), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
#Remove junk_genes
michelle_raw <- michelle_raw[which(rownames(michelle_raw) %in% junk_genes == FALSE),]
#Recalculate p_adj_score
michelle_raw$p_adj_score <- p.adjust(michelle_raw$p_score, "fdr")

#Get names of sbt significant genes
michelle_sig <- rownames(michelle_raw[which(michelle_raw$p_adj_score <= 0.05),])
my_sig <- rownames(sbt_raw[which(sbt_raw$adj_p_score <= 0.05),])
#Get overlap
overlap <- intersect(my_sig, michelle_sig)
#The overlap is 100% of my significant results and 72% of her set

#Clean-up junk
remove(michelle_sig, michelle_raw)

#Calculate intersection with overall_analysis set
overlap <- intersect(rownames(overall_raw[which(overall_raw$adj_p_score < 0.05),]), my_sig)

# 3) Read in files for X-linked-GRM df change ====

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

# 4) Read in files for autosome-GRM df change ====

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
remove(covar, i)

#Remove junk genes
junk_genes <- c("HLA-DPA1", "HLA-DRB5", "HLA-G")
foldchange_df_auto <- foldchange_df_auto[rownames(foldchange_df_auto) %in% junk_genes == FALSE,]
#Do same for expr_final but also remove non-HLA genes
expr_final <- expr_final[which(rownames(expr_final) %in% junk_genes == FALSE & substr(rownames(expr_final), 1, 4) == "HLA-"),]

#Order of these next two lines is important!!!!!
#Merge x-GRM for x-linked genes with auto-GRM for autosomal genes into one data frame
foldchange_df_x <- rbind(cbind(foldchange_df_x, "Chr_type"=rep("X-linked", nrow(foldchange_df_x))), cbind(foldchange_df_auto, "Chr_type"=rep("Autosomal", nrow(foldchange_df_auto))))
#Merge two autosomal-GRM-based foldchange_dfs into one df
foldchange_df <- rbind(cbind(foldchange_df, "Chr_type"=rep("X-linked", nrow(foldchange_df))), cbind(foldchange_df_auto, "Chr_type"=rep("Autosomal", nrow(foldchange_df_auto))))


#Remove junk data
suppressWarnings(remove(auto_raw, pheno_HLA, HLA_individuals_bound, combined, export_list, individuals_bound, pvals, foldchange_df_x, foldchange_df_auto, auto_raw_auto))



# 5) Make overall foldchange data ====

#Set directories
inp_dir <- "C:/Users/Janel/Documents/sex-specific/processed_data/"
filter_dir <- "C:/Users/Janel/Documents/sex-specific/filters/"

#Filter foldchange_df just for autosomal genes since that's all we need for the figures
foldchange_df_autosomal <- foldchange_df[which(foldchange_df$Chr_type == "Autosomal"),]
#Sort dataframes to prep for merging
foldchange_df_autosomal <- foldchange_df_autosomal[order(rownames(foldchange_df_autosomal)),]
overall_raw <- overall_raw[order(rownames(overall_raw)),]
sbt_raw <- sbt_raw[order(rownames(sbt_raw)),]
#Bind adj_p_score columns to foldchange_df_autosomal
foldchange_df_autosomal <- cbind.data.frame(foldchange_df_autosomal, p_adj_overall=overall_raw$adj_p_score, p_adj_sbt=sbt_raw$adj_p_score)


#Read in id files
males_iid_file <- read.table(paste(filter_dir, "males.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
males_HLA_iid_file <- read.table(paste(filter_dir, "HLA.males.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
females_iid_file <- read.table(paste(filter_dir, "females.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
females_HLA_iid_file <- read.table(paste(filter_dir, "HLA.females.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
iid_file <- rbind.data.frame(males_iid_file, females_iid_file)
HLA_iid_file <- rbind.data.frame(males_HLA_iid_file, females_HLA_iid_file)
remove(females_HLA_iid_file, males_HLA_iid_file)

#Read in list of HLA genes
HLA_genes <- read.table(paste(filter_dir, "HLA_auto_genes_DE.tsv", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
HLA_genes <- as.character(HLA_genes$V1)

#Read in expression
expr_raw <- read.table(paste(inp_dir, "BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#sort by row names
expr_raw <- expr_raw[order(row.names(expr_raw)),]
#Exclude junk genes
expr_raw <- expr_raw[which(rownames(expr_raw) %in% junk_genes == FALSE),]
#Sort expr_raw
expr_raw <- expr_raw[order(rownames(expr_raw)),]

#Read in pheno files
pheno_males <- read.table(paste(inp_dir, "males.Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pheno_females <- read.table(paste(inp_dir, "females.Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pheno_raw <- rbind.data.frame(pheno_females[,1:20], pheno_males[,1:20])
remove(pheno_males, pheno_females)

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
foldchange_df_autosomal<-cbind(foldchange_df_autosomal, foldchange_auto)

#Clear junk
suppressWarnings(remove(expr_raw, pheno_HLA, pheno_raw, foldchange_auto, gender, genders))

# 6) Run fold change tests ====

#Create foldchange overlap dfs
foldchange_df_overall <- foldchange_df[rownames(overall_raw[which(overall_raw$adj_p_score <= 0.05),]),]
foldchange_df_interaction <- foldchange_df[rownames(sbt_raw[which(sbt_raw$adj_p_score <= 0.05),]),]

#Append abs foldchange columns
foldchange_df_overall <- cbind.data.frame(foldchange_df_overall, females_FC_abs=abs(foldchange_df_overall$females_FC), males_FC_abs=abs(foldchange_df_overall$males_FC))
foldchange_df_interaction <- cbind.data.frame(foldchange_df_interaction, females_FC_abs=abs(foldchange_df_interaction$females_FC), males_FC_abs=abs(foldchange_df_interaction$males_FC))

#Run Wilcoxon tests
overall_test <- wilcox.test(foldchange_df_overall[,"females_FC_abs"], foldchange_df_overall[,"males_FC_abs"], paired = TRUE, alternative = "two.sided")
interaction_test <- wilcox.test(foldchange_df_interaction[,"females_FC_abs"], foldchange_df_interaction[,"males_FC_abs"], paired = TRUE, alternative = "two.sided")

#Calculate numbers and percents of genes where males_FC_abs > females_FC_abs
num_interaction <- foldchange_df_interaction[which(foldchange_df_interaction$males_FC_abs > foldchange_df_interaction$females_FC_abs),]
num_overall <- foldchange_df_overall[which(foldchange_df_overall$males_FC_abs > foldchange_df_overall$females_FC_abs),]
pct_interaction <- nrow(num_interaction)/nrow(foldchange_df_interaction)
pct_overall <- nrow(num_overall)/nrow(foldchange_df_overall)


# 7) Make plots for manuscript ====

#clean junk
suppressWarnings(remove(num_overall, num_interaction, foldchange_df_interaction, foldchange_df_overall, gender, genders, num_sbt, overlap))

#cbind auto_sig
foldchange_df_autosomal <- cbind.data.frame(foldchange_df_autosomal, 
                                            neg_log_overall=-log(foldchange_df_autosomal$p_adj_overall),
                                            overall_sig=ifelse(foldchange_df_autosomal$p_adj_overall <= 0.05, "Sig", "Not Sig"))

#Set out directory
out_dir <- "C:/Users/Janel/Documents/sex-specific/overall_analyses/"
#Make volcano plot
pdf(paste(out_dir, "Supplemental_Figure_1.pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
suppressWarnings(print(ggplot(foldchange_df_autosomal[which(foldchange_df_autosomal$Chr_type == "Autosomal"),], aes(x = foldchange_auto, y=neg_log_overall)) + #volcanoplot with log2Foldchange versus pvalue
                         geom_point(aes(col=overall_sig)) + #add points colored by significance
                         scale_color_manual(values=c("black", "red")) + 
                         ylim(0,18) + 
                         xlim(-1.5,1.5) + 
                         theme_bw() + 
                         theme(axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank()) + 
                         theme(legend.position="none", legend.title = element_blank(), plot.title = element_blank()) + 
                         xlab("log(fold change)") + ylab(expression(paste(-log[10],"(pvalue)",sep="")))))
dev.off()

#Append color column to foldchange_df_autosomal
foldchange_df_autosomal <- cbind.data.frame(foldchange_df_autosomal, 
                                            scatter_color=ifelse(abs(foldchange_df_autosomal$males_FC) > abs(foldchange_df_autosomal$females_FC), "blue", "red"))
#Append name column to foldchange_df_autosomal
foldchange_df_autosomal <- cbind.data.frame(foldchange_df_autosomal, 
                                            scatter_label=ifelse(foldchange_df_autosomal$males_FC > 6 |
                                                                   foldchange_df_autosomal$females_FC > 6 |
                                                                   foldchange_df_autosomal$males_FC < -4 |
                                                                   foldchange_df_autosomal$females_FC < -4, rownames(foldchange_df_autosomal), ""))

#Make scatter plot
pdf(paste(out_dir, "Figure_4.pdf", sep = ""), width = 8, height = 8, useDingbats = FALSE)
suppressWarnings(ggplot(foldchange_df_autosomal[which(foldchange_df_autosomal$Chr_type == "Autosomal" 
                                                      & rownames(foldchange_df_autosomal) %in% my_sig),], #Filter for 1,055 sbt genes
                        aes(x = females_FC, y = males_FC)) + 
                   geom_point(aes(col=scatter_color)) + #add points colored by significance and nobody cares about escape genes
                   scale_color_manual(values=c("cornflowerblue", "coral2")) +
                   ylim(-8,10.5) + 
                   xlim(-8,10.5) + 
                   #theme_classic() + 
                   theme(legend.position="none", legend.title = element_blank()) +
                   geom_text_repel(aes(label=scatter_label), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 0.005) + 
                   geom_hline(yintercept = 0, linetype="solid", color="black") + 
                   geom_vline(xintercept = 0, linetype="solid", color="black") + 
                   geom_abline(slope = 1, intercept = 0, size=0.1, color="black", alpha = 0.2) + 
                   xlab("Female response to LPS \n (fold change)") + ylab("Male response to LPS \n (fold change)"))
dev.off()


# 8) Check HLA Enrichment ====

#Get gene list names
expressed_autosomal <- rownames(overall_raw)
expressed_HLA <- expressed_autosomal[which(regexpr("HLA", expressed_autosomal) > 0 & regexpr("HHLA", expressed_autosomal) == -1)]

#Make sbt enrichment table and test
sbt_table <- cbind.data.frame(HLA=c(length(expressed_HLA) - length(intersect(expressed_HLA, my_sig)), length(intersect(expressed_HLA, my_sig))), 
                              ALL=c(length(expressed_autosomal) - length(intersect(expressed_autosomal, my_sig)), length(intersect(expressed_autosomal, my_sig))))
sbt_test <- fisher.test(sbt_table)

#Make overall table and test
overall_table <- cbind.data.frame(HLA=c(length(expressed_HLA) - nrow(overall_raw[which(rownames(overall_raw) %in% expressed_HLA & overall_raw$adj_p_score < 0.05),]), 
                                        nrow(overall_raw[which(rownames(overall_raw) %in% expressed_HLA & overall_raw$adj_p_score < 0.05),])), 
                                  ALL=c(length(expressed_autosomal) - nrow(overall_raw[which(rownames(overall_raw) %in% expressed_autosomal & overall_raw$adj_p_score < 0.05),]), 
                                        nrow(overall_raw[which(rownames(overall_raw) %in% expressed_autosomal & overall_raw$adj_p_score < 0.05),])))
overall_test <- fisher.test(overall_table)

#Write HLA file
HLA_output <- foldchange_df_autosomal[expressed_HLA, c("foldchange_auto", "females_FC", "males_FC", "p_adj_overall", "p_adj_sbt")]
write.table(HLA_output, paste(out_dir, "HLA_table_1.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

# 9) Make lists for interferon enrichment test ====

#Make gene lists
sbt_genes <- rownames(sbt_raw[which(sbt_raw$adj_p_score < 0.05),])
not_sbt_genes <- rownames(overall_raw[rownames(overall_raw) %in% rownames(sbt_raw[which(sbt_raw$adj_p_score >= 0.05),]),])
overall_but_not_sbt_genes <- rownames(overall_raw[which(overall_raw$adj_p_score < 0.05 & rownames(overall_raw) %in% sbt_genes == FALSE),])
#Write to file
write.table(sbt_genes, file = paste(out_dir, "sbt_genes.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(not_sbt_genes, file = paste(out_dir, "not_sbt_genes.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(overall_but_not_sbt_genes, file = paste(out_dir, "overall_but_not_sbt_genes.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#Clear junk
remove(sbt_genes, not_sbt_genes, overall_but_not_sbt_genes)

# 10) Make HLA boxplots ====

#Clean junk
remove(foldchange_df, expressed_autosomal, pct_interaction, pct_overall, sbt_test, sbt_table, overall_test, overall_table)

#Bind sex info to the two individuals_bound table
HLA_individuals_bound <- cbind.data.frame(HLA_individuals_bound, sex=ifelse(HLA_individuals_bound[,1] %in% females_iid_file[,2], "Female", "Male"))
individuals_bound <- cbind.data.frame(individuals_bound, sex=ifelse(individuals_bound[,1] %in% females_iid_file[,2], "Female", "Male"))
#Make sex_tables for function to use
HLA_sex_table <- rbind.data.frame(cbind.data.frame(extract_ord=HLA_individuals_bound[,2], sex=HLA_individuals_bound[,4]),
                                  cbind.data.frame(extract_ord=HLA_individuals_bound[,3], sex=HLA_individuals_bound[,4]))
sex_table <- rbind.data.frame(cbind.data.frame(extract_ord=individuals_bound[,2], sex=individuals_bound[,4]),
                                  cbind.data.frame(extract_ord=individuals_bound[,3], sex=individuals_bound[,4]))
remove(HLA_individuals_bound, individuals_bound)
rownames(HLA_sex_table) <- HLA_sex_table[,1]
rownames(sex_table) <- sex_table[,1]

#Read in grant iids to append treatment info to sex table
untreated_iids <- read.table(paste(filter_dir, "grant.iids.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
untreated_iids <- untreated_iids[,2]
sex_table <- cbind.data.frame(sex_table, treat=ifelse(as.numeric(rownames(sex_table)) %in% untreated_iids, "VEH", "LPS"))
suppressWarnings(remove(untreated_iids, females_iid_file, males_iid_file, templist))

#Define out_directory
out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/HLA_boxplots/"

#Define function to make boxplot given 
make_box_plot <- function(sex_table, expr_final, out_dir, gene, p_value, fig_number) {
  
  #transpose expr_final and bind on treatment and sex info
  t_expr_final <- as.data.frame(t(expr_final))
  t_expr_final <- cbind.data.frame(t_expr_final, 
                                   sex=sex_table[substr(rownames(t_expr_final), 2, length(rownames(t_expr_final))), 2], 
                                   treat=sex_table[substr(rownames(t_expr_final), 2, length(rownames(t_expr_final))), 3])
  t_expr_final <- cbind.data.frame(t_expr_final, combo=paste(t_expr_final$sex, "/", t_expr_final$treat, sep = ""))
  #bind on color column
  t_expr_final <- cbind.data.frame(t_expr_final, color=ifelse(t_expr_final$combo == "Female/LPS", "red", 
                                                              ifelse(t_expr_final$combo == "Female/VEH", "coral2", 
                                                                     ifelse(t_expr_final$combo == "Male/LPS", "blue", "skyblue"))))
  

  #Make boxplots
  #open jpg
  jpeg(paste(out_dir, gene, ".jpg", sep = ""))
  #make graphs
  print(ggplot(data=t_expr_final, aes(x=combo, y=get(gene))) + geom_boxplot(color = c("red", "coral2", "blue", "skyblue"), outlier.shape = NA) +
          labs(x="Sex/Treatment", y=expression(paste(log[2],"(expression)",sep=""))) + 
          ggtitle(paste(gene, " Expression", sep = "")) + 
          theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold.italic")) +
          geom_jitter(shape=16, position=position_jitter(0.1), color = t_expr_final$color, show.legend = FALSE) +  
          annotate("text", -Inf, Inf, label = paste("FDR adj.p (interaction test) = ", formatC(p_value, format = "f", digits = 3), sep = ""), hjust = -0.05, vjust = 1, size = 6))
          
  #close jpg
  dev.off()
  
  #Make pdf boxplots
  if(gene == "HLA-DQA1" | gene == "HLA-DQB1"){
    #open pdf
    pdf(paste(out_dir, "Figure_", fig_number, ".pdf", sep = ""))
    #make graphs
    print(ggplot(data=t_expr_final, aes(x=combo, y=get(gene))) + geom_boxplot(color = c("red", "coral2", "blue", "skyblue"), outlier.shape = NA) +
            labs(x="Sex/Treatment", y=expression(paste(log[2],"(expression)",sep=""))) +
            theme(plot.title = element_blank()) +
            geom_jitter(shape=16, position=position_jitter(0.1), color = t_expr_final$color, show.legend = FALSE) +  
            annotate("text", -Inf, Inf, label = paste("FDR adj.p (interaction test) = ", formatC(p_value, format = "f", digits = 3), sep = ""), hjust = -0.05, vjust = 1, size = 6))
    
    #close jpg
    dev.off()
  } else {
    #open pdf
    pdf(paste(out_dir, "Figure_", fig_number, ".pdf", sep = ""))
    #make graphs
    print(ggplot(data=t_expr_final, aes(x=combo, y=get(gene))) + geom_boxplot(color = c("red", "coral2", "blue", "skyblue"), outlier.shape = NA) +
            labs(x="Sex/Treatment", y=expression(paste(log[2],"(expression)",sep=""))) +
            theme(plot.title = element_blank()) +
            geom_jitter(shape=16, position=position_jitter(0.1), color = t_expr_final$color, show.legend = FALSE) +  
            annotate("text", -Inf, Inf, label = paste("FDR adj.p (interaction test) = ", formatC(p_value, format = "f", digits = 3), sep = ""), hjust = -0.05, vjust = 1, size = 6))
    
    #close jpg
    dev.off()
  }
  
}

#Call function
make_box_plot(sex_table, expr_final, out_dir, "HLA-A", sbt_raw["HLA-A", "adj_p_score"], "7A")
make_box_plot(sex_table, expr_final, out_dir, "HLA-B", sbt_raw["HLA-B", "adj_p_score"], "7B")
make_box_plot(sex_table, expr_final, out_dir, "HLA-DQA1", sbt_raw["HLA-DQA1", "adj_p_score"], "7C")
make_box_plot(sex_table, expr_final, out_dir, "HLA-DQB1", sbt_raw["HLA-DQB1", "adj_p_score"], "7D")





# 11) Check Piasecka Genes ====

#Read in gene lists
Piasecka_all_genes <- read.table("C:/Users/Janel/Documents/sex-specific/overall_analyses/Piasecka_all_genes.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Piasecka_all_genes <- Piasecka_all_genes[,1]

#Calculate sbt genes in overlap
length(intersect(Piasecka_all_genes, my_sig))

