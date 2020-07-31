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
library(onewaytests)


# 1) Function MakeCombinedDataFrame(females_auto, males_auto, females_x, males_x) ====

MakeCombinedDataFrame <- function(females_auto, males_auto, females_x, males_x){
  
  #Make a single variance data frame for each gender
  females <- rbind(cbind(females_auto, rep("auto", length(females_auto))), cbind(females_x, rep("x", length(females_x))))
  males <- rbind(cbind(males_auto, rep("auto", length(males_auto))), cbind(males_x, rep("x", length(males_x))))
  colnames(females) <- c("females_stat", "chr_type")
  colnames(males) <- c("males_stat", "chr_type")
  #Bind Gender columns
  females <- cbind.data.frame(females, gender=rep("females", nrow(females)))
  males <- cbind.data.frame(males, gender=rep("males", nrow(males)))
  #Remove non-overlapping rows and sort
  females <- females[which(rownames(females) %in% rownames(males)),]
  males <- males[which(rownames(males) %in% rownames(females)),]
  females <- females[sort(rownames(females)),]
  males <- males[sort(rownames(males)),]
  #Make Gene a column instead of rownames
  females <- cbind.data.frame(gene=as.character(rownames(females)), females)
  males <- cbind.data.frame(gene=as.character(rownames(males)), males)
  #Bind Data frames together
  combined_df <- cbind.data.frame(gene=females$gene, females_stat=females$females_stat, males_stat=males$males_stat, chr_type=females$chr_type)
  #Convert variance column to numeric data
  combined_df$females_stat <- as.numeric(as.character(combined_df$females_stat))
  combined_df$males_stat <- as.numeric(as.character(combined_df$males_stat))
  combined_df$gene <- as.character(combined_df$gene)
  return(combined_df)
  
}

# 2) Read in files for df change ====

#Set input directory
inp_dir <- "C:/Users/Janel/Documents/sex-specific/processed_data/"
filter_dir <- "C:/Users/Janel/Documents/sex-specific/filters/"

#Define genders
genders <- c("females", "males")

#Create empty data frames to be appended to
expr_final <- as.data.frame(matrix(nrow = 393))
covar <- as.data.frame(matrix(ncol = 7))
colnames(covar) <- c("Extract_Ord", "sex", "treat", "Age", "Pool_1.0", "library_tech", "FINDIV")
foldchange_df<-as.data.frame(matrix(nrow = 393))
foldchange_df_x<-as.data.frame(matrix(nrow = 393))

#Repeat for both genders
for (gender in genders) {
  
  #Read in X_GRM id list
  x_grm_ids <- read.table(paste(filter_dir, gender, ".X.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  #Read in list of ids
  iid_file <- read.table(paste(filter_dir, gender, ".iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  #Read in expression
  expr_raw <- read.table(paste(inp_dir, gender, ".logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  #sort by row names
  expr_raw <- expr_raw[order(row.names(expr_raw)),]
  
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
  #Filter pheno file
  pheno_raw <- pheno_raw[which(pheno_raw$FINDIV %in% iid_file$V2),]
  
  #Get list of ids
  individuals<-pheno_raw[,"FINDIV"]
  individuals<-subset(individuals,duplicated(individuals)!=TRUE)
  individuals_bound<-cbind(individuals,individuals,individuals)
  for (k in 1:nrow(individuals_bound)) {
    individuals_bound[k,2]<-pheno_raw[which(pheno_raw$FINDIV==individuals_bound[k,1] & pheno_raw$treat=="LPS"),][,"Extract_Ord"]
    individuals_bound[k,3]<-pheno_raw[which(pheno_raw$FINDIV==individuals_bound[k,1] & pheno_raw$treat=="null"),][,"Extract_Ord"]
  }
  remove(individuals,k,iid_file)
  colnames(individuals_bound) <- c("FINDIV", "LPS", "null")
  
  #Filter individuals_bound for the X_GRM
  individuals_bound_x <- individuals_bound[which(individuals_bound[,2] %in% x_grm_ids[,2]),]
  
  #Define fold change column
  foldchange<-vector()
  foldchange_x<-vector()
  #auto_grm
  for (i in 1:nrow(expr_raw)) {
    
    #Create fold change column (method 1 divide then average)
    templist<-vector()
    for (j in 1:nrow(individuals_bound)) {
      templist[j]<-expr_raw[i,sprintf("X%04d",individuals_bound[j,2])] - expr_raw[i,sprintf("X%04d",individuals_bound[j,3])]
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
      templist[j]<-expr_raw[i,sprintf("X%04d",individuals_bound_x[j,2])] - expr_raw[i,sprintf("X%04d",individuals_bound_x[j,3])]
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
  
  #Store individuals bound
  if (gender == "females") {
    females_key <- individuals_bound 
  } else if (gender == "males") {
    males_key <- individuals_bound 
  }
  
  
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
suppressWarnings(remove(pheno_raw, expr_raw,templist, foldchange, foldchange_x, individuals_bound, individuals_bound_x, x_grm_ids))

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

# 3) Read in expression files again ====

#Remove foldchange data
remove(foldchange_df, foldchange_df_x, gender, genders)

#Create list of junk genes
junk_genes <- c("HLA-DPA1", "HLA-DRB5", "HLA-G")

#Read in HLA iids
females_HLA_iids <- read.table(paste(filter_dir, "HLA.females.iids.txt", sep = ""), sep = "\t", stringsAsFactors = FALSE)
males_HLA_iids <- read.table(paste(filter_dir, "HLA.males.iids.txt", sep = ""), sep = "\t", stringsAsFactors = FALSE)

#Make HLA_keys
females_HLA_key <- females_key[which(females_key[,1] %in% females_HLA_iids$V2),]
males_HLA_key <- males_key[which(males_key[,1] %in% males_HLA_iids$V2),]
remove(females_HLA_iids, males_HLA_iids)

#Read in x expr files
expr_males_x <- read.table(paste(inp_dir, "males", ".logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
expr_females_x <- read.table(paste(inp_dir, "females", ".logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#Filter for ids in keys
expr_females_x <- cbind.data.frame(expr_females_x[,sprintf("X%04d", females_key[,2])], expr_females_x[,sprintf("X%04d", females_key[,3])])
expr_males_x <- cbind.data.frame(expr_males_x[,sprintf("X%04d", males_key[,2])], expr_males_x[,sprintf("X%04d", males_key[,3])])

#Read in auto expr files
expr_auto <- read.table(paste(inp_dir, "BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#Filter out junk genes
expr_auto <- expr_auto[which(rownames(expr_auto) %in% junk_genes == FALSE),]
#Split expression by sex
expr_males_auto <- cbind.data.frame(expr_auto[,paste("X", males_key[,2], sep = "")], expr_auto[,paste("X", males_key[,3], sep = "")])
expr_females_auto <- cbind.data.frame(expr_auto[,paste("X", females_key[,2], sep = "")], expr_auto[,paste("X", females_key[,3], sep = "")])
remove(expr_auto)

#Get data set overlaps
expr_males_x <- expr_males_x[which(rownames(expr_males_x) %in% rownames(expr_females_x)),]
expr_females_x <- expr_females_x[which(rownames(expr_females_x) %in% rownames(expr_males_x)),]
expr_males_auto <- expr_males_auto[which(rownames(expr_males_auto) %in% rownames(expr_females_auto)),]
expr_females_auto <- expr_females_auto[which(rownames(expr_females_auto) %in% rownames(expr_males_auto)),]

#Overwrite colnames of expr_x matrices 
colnames(expr_females_x) <- colnames(expr_females_auto)
colnames(expr_males_x) <- colnames(expr_males_auto)

#Read in HLA data to overwrite rows of auto expression matrices
expr_HLA <- read.table(paste(inp_dir, "HLA_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#Separate into two data frames by sex
expr_females_HLA <- cbind.data.frame(expr_HLA[,paste("X", females_HLA_key[,2], sep = "")], expr_HLA[,paste("X", females_HLA_key[,3], sep = "")])
expr_males_HLA <- cbind.data.frame(expr_HLA[,paste("X", males_HLA_key[,2], sep = "")], expr_HLA[,paste("X", males_HLA_key[,3], sep = "")])
remove(expr_HLA)


# 4) Function ExtractExpr(gene_expr_data, sex_key) ====

ExtractExpr <- function(gene_expr_data, sex_key, treat){
  
  #Calculate Gene Response
  expression_of_type <- mean(gene_expr_data[paste("X",sex_key[,treat], sep = "")])
  #Return value
  return(expression_of_type)
  
}

# 5) Compare Expression Levels for Emma (Commented Out) ====

##Make two data frames, one for LPS and one for VEH expression
##We'll need to compare these to understand whether males have lower starting expression levels
#avg_expr_LPS <- MakeCombinedDataFrame(apply(expr_females_auto, 1, ExtractExpr, sex_key=females_key, treat="LPS"),
#                                      apply(expr_males_auto, 1, ExtractExpr, sex_key=males_key, treat="LPS"),
#                                      apply(expr_females_x, 1, ExtractExpr, sex_key=females_key, treat="LPS"),
#                                      apply(expr_males_x, 1, ExtractExpr, sex_key=males_key, treat="LPS"))
#avg_expr_VEH <- MakeCombinedDataFrame(apply(expr_females_auto, 1, ExtractExpr, sex_key=females_key, treat="null"),
#                                      apply(expr_males_auto, 1, ExtractExpr, sex_key=males_key, treat="null"),
#                                      apply(expr_females_x, 1, ExtractExpr, sex_key=females_key, treat="null"),
#                                      apply(expr_males_x, 1, ExtractExpr, sex_key=males_key, treat="null"))
#
##Reorganize dfs by chromosome type
#avg_expr_x <- cbind.data.frame(gene=avg_expr_LPS$gene, chr_type=avg_expr_LPS$chr_type, 
#                             females_VEH=avg_expr_VEH$females_stat, females_LPS=avg_expr_LPS$females_stat,
#                             males_VEH=avg_expr_VEH$males_stat, males_LPS=avg_expr_LPS$males_stat)
#avg_expr_auto <- avg_expr_x[which(avg_expr_x$chr_type == "auto"),]
#avg_expr_x <- avg_expr_x[which(avg_expr_x$chr_type == "x"),]
##Convert factor columns to characters
#avg_expr_auto$gene <- as.character(avg_expr_auto$gene)
#avg_expr_x$gene <- as.character(avg_expr_x$gene)

##Add columns showing males and females FC
#avg_expr_x$females_FC <- avg_expr_x$females_LPS - avg_expr_x$females_VEH
#avg_expr_x$males_FC <- avg_expr_x$males_LPS - avg_expr_x$males_VEH
#avg_expr_auto$females_FC <- avg_expr_auto$females_LPS - avg_expr_auto$females_VEH
#avg_expr_auto$males_FC <- avg_expr_auto$males_LPS - avg_expr_auto$males_VEH


## prep data class (goes up or down)
#avg_expr_x$females_class <- ifelse((avg_expr_x$females_VEH - avg_expr_x$females_LPS) >= 0, "#E95C20FF", "green")
#avg_expr_x$males_class <- ifelse((avg_expr_x$males_VEH - avg_expr_x$males_LPS) >= 0, "#E95C20FF", "green")
#avg_expr_auto$females_class <- ifelse((avg_expr_auto$females_VEH - avg_expr_auto$females_LPS) >= 0, "#E95C20FF", "green")
#avg_expr_auto$males_class <- ifelse((avg_expr_auto$males_VEH - avg_expr_auto$males_LPS) >= 0, "#E95C20FF", "green")

## Recombine data frames to work for box plot
#auto_boxplot_df <- cbind.data.frame(gene=as.character(append(avg_expr_auto$gene, append(avg_expr_auto$gene,
#                                                     append(avg_expr_auto$gene, avg_expr_auto$gene)))),
#                                    expression=append(avg_expr_auto$females_VEH, append(avg_expr_auto$females_LPS, 
#                                                     append(avg_expr_auto$males_VEH, avg_expr_auto$males_LPS))),
#                                    category=append(rep("Females/VEH", nrow(avg_expr_auto)), append(rep("Females/LPS", nrow(avg_expr_auto)),
#                                                   append(rep("Males/VEH", nrow(avg_expr_auto)), rep("Males/LPS", nrow(avg_expr_auto))))),
#                                    FC_color=as.character(append(avg_expr_auto$females_class, append(avg_expr_auto$females_class, 
#                                                   append(avg_expr_auto$males_class, avg_expr_auto$males_class)))))
#x_boxplot_df <- cbind.data.frame(gene=as.character(append(avg_expr_x$gene, append(avg_expr_x$gene,
#                                                                                        append(avg_expr_x$gene, avg_expr_x$gene)))),
#                                    expression=append(avg_expr_x$females_VEH, append(avg_expr_x$females_LPS, 
#                                                                                        append(avg_expr_x$males_VEH, avg_expr_x$males_LPS))),
#                                    category=append(rep("Females/VEH", nrow(avg_expr_x)), append(rep("Females/LPS", nrow(avg_expr_x)),
#                                                                                                    append(rep("Males/VEH", nrow(avg_expr_x)), rep("Males/LPS", nrow(avg_expr_x))))),
#                                    FC_color=as.character(append(avg_expr_x$females_class, append(avg_expr_x$females_class, 
#                                                                                                     append(avg_expr_x$males_class, avg_expr_x$males_class)))))
##Set out_dir
#out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/"
## Plot auto
#tiff(paste(out_dir, "autosomal_boxplot.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
#print(ggplot(data=auto_boxplot_df, aes(x=category, y=expression)) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
#        labs(x="", y="Log-Normalized Expression") + 
#        theme(plot.title = element_blank()) + 
#        geom_jitter(shape=16, position=position_jitter(0.1), show.legend = FALSE, aes(color = FC_color)) +
#        scale_color_manual(values = c("green", "#E95C20FF"))
#      )
#dev.off()
##Plot x
#tiff(paste(out_dir, "x_linked_boxplot.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
#print(ggplot(data=x_boxplot_df, aes(x=category, y=expression)) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
#        labs(x="", y="Log-Normalized Expression") + 
#        theme(plot.title = element_blank()) + 
#        geom_jitter(shape=16, position=position_jitter(0.1), show.legend = FALSE, aes(color = FC_color)) +
#        scale_color_manual(values = c("green", "#E95C20FF"))
#)
#dev.off()


# Plot Change Graph (Commented Out because it's too messy) 

#ggplot(avg_expr_x) + geom_segment(aes(x=1, xend=5, y=males_VEH, yend=males_LPS, col=males_class), size=.75, show.legend=F) +
#  geom_vline(xintercept=1, linetype="dashed", size=.1) + 
#  geom_vline(xintercept=5, linetype="dashed", size=.1) +
#  geom_hline(yintercept=0, linetype="dashed", size=.1) +
#  scale_fill_manual(values = c(alpha("gray", alpha = 0.1), "black")) + 
#  scale_color_manual(labels = c("Down", "Up"), 
#                     values = c("#E95C20FF"="#f8766d", "green"="#00ba38")) +  # color of lines and points
#  labs(x="", y="log(Normalized Expression)") +  # Axis labels
#  xlim(1, 5) + ylim(min(min(avg_expr_x$males_LPS),min(avg_expr_x$males_VEH), min(avg_expr_x$males_LPS), min(avg_expr_x$males_VEH)),
#                    max(max(avg_expr_x$males_LPS),max(avg_expr_x$males_VEH), max(avg_expr_x$males_LPS), max(avg_expr_x$males_VEH))) + # X and Y axis limits
#  theme(panel.background = element_blank(), 
#        panel.grid = element_blank(),
#        axis.ticks = element_blank(),
#        axis.text.x = element_blank(),
#        panel.border = element_blank(),
#        plot.margin = unit(c(1,2,1,2), "cm"))

# 5.5) Calculate some basic stats for Emma (Commented Out) ====

##Tables of ups and downs
#table(avg_expr_x$females_class, avg_expr_x$males_class) #female categories are displayed on left side
#table(avg_expr_auto$females_class, avg_expr_auto$males_class) #female categories are displayed on left side
##Average Levels
#average_expression_x <-c(females_VEH=mean(avg_expr_x$females_VEH), females_LPS=mean(avg_expr_x$females_LPS), 
#                         males_VEH=mean(avg_expr_x$males_VEH), males_LPS=mean(avg_expr_x$males_LPS))
#average_expression_auto <-c(females_VEH=mean(avg_expr_auto$females_VEH), females_LPS=mean(avg_expr_auto$females_LPS), 
#                            males_VEH=mean(avg_expr_auto$males_VEH), males_LPS=mean(avg_expr_auto$males_LPS))
##Median Levels
#median_expression_x <-c(females_VEH=median(avg_expr_x$females_VEH), females_LPS=median(avg_expr_x$females_LPS), 
#                        males_VEH=median(avg_expr_x$males_VEH), males_LPS=median(avg_expr_x$males_LPS))
#median_expression_auto <-c(females_VEH=median(avg_expr_auto$females_VEH), females_LPS=median(avg_expr_auto$females_LPS), 
#                           males_VEH=median(avg_expr_auto$males_VEH), males_LPS=median(avg_expr_auto$males_LPS))
##Average Response Changes
#avg_FC_females_x <- mean(avg_expr_x$females_LPS - avg_expr_x$females_VEH)
#avg_FC_males_x <- mean(avg_expr_x$males_LPS - avg_expr_x$males_VEH)
#avg_FC_females_auto <- mean(avg_expr_auto$females_LPS - avg_expr_auto$females_VEH)
#avg_FC_males_auto <- mean(avg_expr_auto$males_LPS - avg_expr_auto$males_VEH)
##Median Response Changes
#med_FC_females_x <- median(avg_expr_x$females_LPS - avg_expr_x$females_VEH)
#med_FC_males_x <- median(avg_expr_x$males_LPS - avg_expr_x$males_VEH)
#med_FC_females_auto <- median(avg_expr_auto$females_LPS - avg_expr_auto$females_VEH)
#med_FC_males_auto <- median(avg_expr_auto$males_LPS - avg_expr_auto$males_VEH)
##Average Absolute Value Response Changes
#avg_abs_FC_females_x <- mean(abs(avg_expr_x$females_LPS - avg_expr_x$females_VEH))
#avg_abs_FC_males_x <- mean(abs(avg_expr_x$males_LPS - avg_expr_x$males_VEH))
#avg_abs_FC_females_auto <- mean(abs(avg_expr_auto$females_LPS - avg_expr_auto$females_VEH))
#avg_abs_FC_males_auto <- mean(abs(avg_expr_auto$males_LPS - avg_expr_auto$males_VEH))
##Median Absolute Value Response Changes
#med_abs_FC_females_x <- median(abs(avg_expr_x$females_LPS - avg_expr_x$females_VEH))
#med_abs_FC_males_x <- median(abs(avg_expr_x$males_LPS - avg_expr_x$males_VEH))
#med_abs_FC_females_auto <- median(abs(avg_expr_auto$females_LPS - avg_expr_auto$females_VEH))
#med_abs_FC_males_auto <- median(abs(avg_expr_auto$males_LPS - avg_expr_auto$males_VEH))
##Get top FCs for each category
#top_10_females_auto <- avg_expr_auto[order(avg_expr_auto$females_FC, decreasing = TRUE),][1:10,]
#top_10_females_x <- avg_expr_x[order(avg_expr_x$females_FC, decreasing = TRUE),][1:10,]
#top_10_males_auto <- avg_expr_auto[order(avg_expr_auto$males_FC, decreasing = TRUE),][1:10,]
#top_10_males_x <- avg_expr_x[order(avg_expr_x$males_FC, decreasing = TRUE),][1:10,]
##Get bottom FCs for each category
#bottom_10_females_auto <- avg_expr_auto[order(avg_expr_auto$females_FC, decreasing = FALSE),][1:10,]
#bottom_10_females_x <- avg_expr_x[order(avg_expr_x$females_FC, decreasing = FALSE),][1:10,]
#bottom_10_males_auto <- avg_expr_auto[order(avg_expr_auto$males_FC, decreasing = FALSE),][1:10,]
#bottom_10_males_x <- avg_expr_x[order(avg_expr_x$males_FC, decreasing = FALSE),][1:10,]

# 6) Function MakeVariancePlots(females_auto_var, males_auto_var, females_x_var, males_x_var, escape_genes, treat) ====

MakeVariancePlots <- function(females_auto_var, males_auto_var, females_x_var, males_x_var, escape_genes, treat){

  
  #Make a single variance data frame for each gender
  females_var <- rbind(cbind(females_auto_var, rep("auto", length(females_auto_var))), cbind(females_x_var, rep("x", length(females_x_var))))
  males_var <- rbind(cbind(males_auto_var, rep("auto", length(males_auto_var))), cbind(males_x_var, rep("x", length(males_x_var))))
  colnames(females_var) <- c("females_variance", "chr_type")
  colnames(males_var) <- c("males_variance", "chr_type")
  #Bind Gender columns
  females_var <- cbind.data.frame(females_var, gender=rep("females", nrow(females_var)))
  males_var <- cbind.data.frame(males_var, gender=rep("males", nrow(males_var)))
  #Remove non-overlapping rows and sort
  females_var <- females_var[which(rownames(females_var) %in% rownames(males_var)),]
  males_var <- males_var[which(rownames(males_var) %in% rownames(females_var)),]
  females_var <- females_var[sort(rownames(females_var)),]
  males_var <- males_var[sort(rownames(males_var)),]
  #Make Gene a column instead of rownames
  females_var <- cbind.data.frame(gene=as.character(rownames(females_var)), females_var)
  males_var <- cbind.data.frame(gene=as.character(rownames(males_var)), males_var)
  #Bind Data frames together
  variance_df <- cbind.data.frame(gene=females_var$gene, females_variance=females_var$females_variance, males_variance=males_var$males_variance, chr_type=females_var$chr_type)
  #Convert variance column to numeric data
  variance_df$females_variance <- as.numeric(as.character(variance_df$females_variance))
  variance_df$males_variance <- as.numeric(as.character(variance_df$males_variance))
  variance_df$gene <- as.character(variance_df$gene)
  #Assign rownames as gene names
  rownames(variance_df) <- variance_df$gene
  #Remove junk
  remove(females_var, males_var, females_auto_var, males_auto_var, females_x_var, males_x_var)
  
  #Calculate regression lines for plot
  auto_lm <- lm(variance_df[which(variance_df$chr_type == "auto"),"males_variance"] ~ variance_df[which(variance_df$chr_type == "auto"),"females_variance"] - 1)
  x_lm <- lm(variance_df[which(variance_df$chr_type == "x"),"males_variance"] ~ variance_df[which(variance_df$chr_type == "x"),"females_variance"] - 1)
  
  #Make plot without regression lines
  out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/variance_plots/"
  tiff(paste(out_dir, "DE_Variance_Comparison_", treat, ".jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  print(suppressWarnings(ggplot(variance_df, aes(x = females_variance, y = males_variance)) + 
                     geom_point(aes(col=chr_type)) + #add points colored by significance and nobody cares about escape genes
                     scale_color_manual(values=c(alpha("#E95C20FF", 0.4), "#00ff00")) +
                     ylim(min(variance_df$males_variance),max(variance_df$males_variance)) + 
                     xlim(min(variance_df$females_variance),max(variance_df$females_variance)) + 
                     #theme_classic() + 
                     theme(legend.position="none", legend.title = element_blank()) +
                     #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                     geom_hline(yintercept = 0, linetype="solid", color="black") + 
                     geom_vline(xintercept = 0, linetype="solid", color="black") + 
                     geom_abline(slope = 1, intercept = 0, size=0.1, linetype="dashed") + 
                     #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                     xlab("Female Fold Change Variance") + ylab("Male Fold Change Variance")))
  dev.off()
  #Make plot with regression lines
  tiff(paste(out_dir, "DE_Variance_Comparison_w_lines_", treat, ".jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  print(suppressWarnings(ggplot(variance_df, aes(x = females_variance, y = males_variance)) + 
                     geom_point(aes(col=chr_type)) + #add points colored by significance and nobody cares about escape genes
                     scale_color_manual(values=c(alpha("#E95C20FF", 0.4), "#00ff00")) +
                     ylim(min(variance_df$males_variance),max(variance_df$males_variance)) + 
                     xlim(min(variance_df$females_variance),max(variance_df$females_variance)) + 
                     #theme_classic() + 
                     theme(legend.position="none", legend.title = element_blank()) +
                     #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                     geom_hline(yintercept = 0, linetype="solid", color="black") + 
                     geom_vline(xintercept = 0, linetype="solid", color="black") + 
                     geom_abline(slope = 1, intercept = 0, size=0.1, linetype="dashed") + 
                     geom_abline(slope = x_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color="#00ff00") + #auto regression line
                     geom_abline(slope = auto_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("#E95C20FF", 1)) + #x regression line
                     #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                     xlab("Female Fold Change Variance") + ylab("Male Fold Change Variance")))
  dev.off()
  
  #Create variance_df_escape
  variance_df_escape <- cbind.data.frame(variance_df, escape_status=ifelse(rownames(variance_df) %in% escape_genes, TRUE, FALSE))
  
  #Calculate regression lines for plot
  x_escape_lm <- lm(variance_df_escape[which(variance_df_escape$chr_type == "x" & variance_df_escape$escape_status == TRUE),"males_variance"] ~ variance_df_escape[which(variance_df_escape$chr_type == "x" & variance_df_escape$escape_status == TRUE),"females_variance"] - 1)
  x_inactivated_lm <- lm(variance_df_escape[which(variance_df_escape$chr_type == "x" & variance_df_escape$escape_status == FALSE),"males_variance"] ~ variance_df_escape[which(variance_df_escape$chr_type == "x" & variance_df_escape$escape_status == FALSE),"females_variance"] - 1)
  
  #filter variance_df_escape into separate dfs for graphs
  variance_df_inactivated <- variance_df_escape[(variance_df_escape$chr_type == "auto" | (variance_df_escape$chr_type == "x" & variance_df_escape$escape_status == FALSE)),]
  variance_df_escape <- variance_df_escape[(variance_df_escape$chr_type == "auto" | (variance_df_escape$chr_type == "x" & variance_df_escape$escape_status == TRUE)),]
  
  #Make escape plot
  out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/variance_plots/"
  tiff(paste(out_dir, "DE_Variance_Comparison_", treat, ".escape.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  print(suppressWarnings(ggplot(variance_df_escape, aes(x = females_variance, y = males_variance)) + 
                           geom_point(aes(col=chr_type)) + #add points colored by significance and nobody cares about escape genes
                           scale_color_manual(values=c(alpha("#E95C20FF", 0.4), "#00ff00")) +
                           ylim(min(variance_df_escape$males_variance),max(variance_df_escape$males_variance)) + 
                           xlim(min(variance_df_escape$females_variance),max(variance_df_escape$females_variance)) + 
                           #theme_classic() + 
                           theme(legend.position="none", legend.title = element_blank()) +
                           #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                           geom_hline(yintercept = 0, linetype="solid", color="black") + 
                           geom_vline(xintercept = 0, linetype="solid", color="black") + 
                           geom_abline(slope = 1, intercept = 0, size=0.1, linetype="dashed") + 
                           geom_abline(slope = x_escape_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color="#00ff00") + #auto regression line
                           geom_abline(slope = auto_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("#E95C20FF", 1)) + #x regression line
                           #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                           xlab("Female Fold Change Variance") + ylab("Male Fold Change Variance")))
  dev.off()
  
  #Make inactivated plot
  out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/variance_plots/"
  tiff(paste(out_dir, "DE_Variance_Comparison_", treat, ".inactivated.jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
  print(suppressWarnings(ggplot(variance_df_inactivated, aes(x = females_variance, y = males_variance)) + 
                           geom_point(aes(col=chr_type)) + #add points colored by significance and nobody cares about escape genes
                           scale_color_manual(values=c(alpha("#E95C20FF", 0.4), "#00ff00")) +
                           ylim(min(variance_df_inactivated$males_variance),max(variance_df_inactivated$males_variance)) + 
                           xlim(min(variance_df_inactivated$females_variance),max(variance_df_inactivated$females_variance)) + 
                           #theme_classic() + 
                           theme(legend.position="none", legend.title = element_blank()) +
                           #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                           geom_hline(yintercept = 0, linetype="solid", color="black") + 
                           geom_vline(xintercept = 0, linetype="solid", color="black") + 
                           geom_abline(slope = 1, intercept = 0, size=0.1, linetype="dashed") + 
                           geom_abline(slope = x_inactivated_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color="#00ff00") + #auto regression line
                           geom_abline(slope = auto_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("#E95C20FF", 1)) + #x regression line
                           #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                           xlab("Female Fold Change Variance") + ylab("Male Fold Change Variance")))
  dev.off()
  
  #Create bar plot of differences in variances
  #Bind difference column to df
  #variance_df <- cbind.data.frame(variance_df, femalesMinusmales=(variance_df$females_variance-variance_df$males_variance))
  #Plot
  #barplot(variance_df[which(variance_df$chr_type == "auto"), "femalesMinusmales"], main = "Female Var - Male Var")
  
}

# 7) Function CalculateGeneResponseVar(gene_expr_data, sex_key) ====

CalculateGeneResponseVar <- function(gene_expr_data, sex_key){
  
  #Calculate Gene Response
  expression_response <- var(gene_expr_data[paste("X",sex_key[,"LPS"], sep = "")]-gene_expr_data[paste("X",sex_key[,"null"], sep = "")])
  #Return value
  return(expression_response)
  
}

# 8) Function CalculateGeneResponse(gene_expr_data, sex_key) ====

CalculateGeneResponse <- function(gene_expr_data, sex_key){
  
  #Calculate Gene Response
  expression_response <- mean(abs(gene_expr_data[paste("X",sex_key[,"LPS"], sep = "")]-gene_expr_data[paste("X",sex_key[,"null"], sep = "")]))
  #Return value
  return(expression_response)
  
}

# 9) Read in escape genes and Call functions ====

#Read in escape genes
blood_escape <- read.table("C:/Users/Janel/Documents/sex-specific/de_analysis/WHLBLOOD_escape_genes_121619.txt", header = FALSE, stringsAsFactors = FALSE)
blood_escape <- as.character(blood_escape[,1])


#Apply CalculateGeneResponseVar to each matrix
response_females_auto_var <- apply(expr_females_auto, 1, CalculateGeneResponseVar, sex_key=females_key)
response_females_x_var <- apply(expr_females_x, 1, CalculateGeneResponseVar, sex_key=females_key)
response_males_auto_var <- apply(expr_males_auto, 1, CalculateGeneResponseVar, sex_key=males_key)
response_males_x_var <- apply(expr_males_x, 1, CalculateGeneResponseVar, sex_key=males_key)
#Repeat process for HLA genes which will be used to overwrite previous values
response_females_HLA_var <- apply(expr_females_HLA, 1, CalculateGeneResponseVar, sex_key=females_HLA_key)
response_males_HLA_var <- apply(expr_males_HLA, 1, CalculateGeneResponseVar, sex_key=males_HLA_key)
#Overwrite HLA values and remove junk_genes
response_females_auto_var[names(response_females_HLA_var)] <- response_females_HLA_var
response_males_auto_var[names(response_males_HLA_var)] <- response_males_HLA_var
response_females_auto_var <- response_females_auto_var[which(names(response_females_auto_var) %in% junk_genes == FALSE)]
response_males_auto_var <- response_males_auto_var[which(names(response_males_auto_var) %in% junk_genes == FALSE)]
#Remove junk
remove(response_females_HLA_var, response_males_HLA_var)

#Calculate Expression Variances
females_auto_var <- apply(expr_females_auto, 1, var)
males_auto_var <- apply(expr_males_auto, 1, var)
females_x_var <- apply(expr_females_x, 1, var)
males_x_var <- apply(expr_males_x, 1, var)
#Repeat process for HLA genes which will be used to overwrite previous values
females_HLA_var <- apply(expr_females_HLA, 1, var)
males_HLA_var <- apply(expr_males_HLA, 1, var)
#Overwrite HLA values and remove junk_genes
females_auto_var[names(females_HLA_var)] <- females_HLA_var
males_auto_var[names(males_HLA_var)] <- males_HLA_var
females_auto_var <- females_auto_var[which(names(females_auto_var) %in% junk_genes == FALSE)]
males_auto_var <- males_auto_var[which(names(males_auto_var) %in% junk_genes == FALSE)]
#Remove junk
remove(females_HLA_var, males_HLA_var)

#Apply MakeVariancePlots function to response data
MakeVariancePlots(response_females_auto_var, response_males_auto_var, response_females_x_var, response_males_x_var, blood_escape, "response")

# 10) #Run Wilcoxon Tests ====

#Make necessary df
test_df <- rbind.data.frame(cbind.data.frame(females_var=response_females_auto_var, males_var=response_males_auto_var, Chr_type=rep("Autosomal", length(response_females_auto_var))),
                            cbind.data.frame(females_var=response_females_x_var, males_var=response_males_x_var, Chr_type=rep("X-linked", length(response_females_x_var))))
test_df <- cbind.data.frame(test_df, escape_status=ifelse(rownames(test_df) %in% blood_escape, TRUE, FALSE))

#Run tests
all_test <- wilcox.test(test_df[,"females_var"], test_df[,"males_var"], paired = TRUE, alternative = "two.sided")
autosomal_test <- wilcox.test(test_df[which(test_df$Chr_type == "Autosomal"),"females_var"], 
                              test_df[which(test_df$Chr_type == "Autosomal"),"males_var"], paired = TRUE, 
                              alternative = "two.sided")
x_test <- wilcox.test(test_df[which(test_df$Chr_type == "X-linked"),"females_var"], 
                      test_df[which(test_df$Chr_type == "X-linked"),"males_var"], paired = TRUE, 
                      alternative = "two.sided")
escape_test <- wilcox.test(test_df[which(test_df$Chr_type == "X-linked" & test_df$escape_status == TRUE),"females_var"], 
                           test_df[which(test_df$Chr_type == "X-linked" & test_df$escape_status == TRUE),"males_var"], 
                           paired = TRUE, alternative = "two.sided")
inactivated_test <- wilcox.test(test_df[which(test_df$Chr_type == "X-linked" & test_df$escape_status == FALSE),"females_var"], 
                                test_df[which(test_df$Chr_type == "X-linked" & test_df$escape_status == FALSE),"males_var"], 
                                paired = TRUE, alternative = "two.sided")

# 11) Run sub-sampling for sex, make plots, and rerun Wilcoxon test on sampled data ====

#Make holding variables for sampling results (will need to delete first column later)
response_females_auto_var_sampled <- response_females_auto_var
response_females_x_var_sampled <- response_females_x_var

#Loop through 200 times
for (i in 1:200) {
  #Get random ids to pull
  females_sampled <- females_key[sample(1:nrow(females_key), nrow(males_key)),]
  #make combined id list
  females_sampled_combined <- paste("X", append(females_sampled[,"LPS"], females_sampled[,"null"]), sep = "")
  #filter females data frames for sample ids
  expr_females_auto_sampled <- expr_females_auto[,which(colnames(expr_females_auto) %in% females_sampled_combined)]
  expr_females_x_sampled <- expr_females_x[,which(colnames(expr_females_x) %in% females_sampled_combined)]
  
  #Apply CalculateGeneResponseVar to sampled females results
  response_females_auto_var_sampled <- cbind.data.frame(response_females_auto_var_sampled, apply(expr_females_auto_sampled, 1, CalculateGeneResponseVar, sex_key=females_sampled))
  response_females_x_var_sampled <- cbind.data.frame(response_females_x_var_sampled, apply(expr_females_x_sampled, 1, CalculateGeneResponseVar, sex_key=females_sampled))

}

#Remove first columns of sampled variance dfs
response_females_auto_var_sampled <- response_females_auto_var_sampled[,2:ncol(response_females_auto_var_sampled)]
response_females_x_var_sampled <- response_females_x_var_sampled[,2:ncol(response_females_x_var_sampled)]
#Rename columns
colnames(response_females_auto_var_sampled) <- paste("X", 1:200, sep = "")
colnames(response_females_x_var_sampled) <- paste("X", 1:200, sep = "")
#Take average for each row
response_females_auto_final <- apply(response_females_auto_var_sampled, 1, mean)
response_females_x_final <- apply(response_females_x_var_sampled, 1, mean)

#Clear junk
suppressWarnings(remove(i, females_sampled))

#Apply MakeVariancePlots function to sampled response data
MakeVariancePlots(response_females_auto_final, response_males_auto_var, response_females_x_final, response_males_x_var, blood_escape, "sampled")

#Remake necessary df
test_df_sampled <- rbind.data.frame(cbind.data.frame(females_var=response_females_auto_final, males_var=response_males_auto_var, Chr_type=rep("Autosomal", length(response_females_auto_var))),
                            cbind.data.frame(females_var=response_females_x_final, males_var=response_males_x_var, Chr_type=rep("X-linked", length(response_females_x_var))))
#Rerun autosomal Wilcox Test
autosomal_test_sampled <- wilcox.test(test_df_sampled[which(test_df_sampled$Chr_type == "Autosomal"),"females_var"], 
                              test_df_sampled[which(test_df_sampled$Chr_type == "Autosomal"),"males_var"], paired = TRUE, 
                              alternative = "two.sided")

# 12) Calculate and Plot CVs and Fano Factors ====

#Remove junk
suppressWarnings(remove(females_auto_var, males_auto_var, females_x_var, males_x_var, expr_females_x_sampled, expr_females_auto_sampled, response_females_auto_var_sampled, response_females_x_var_sampled))

#Apply CalculateGeneResponse to each matrix
response_females_auto <- apply(expr_females_auto, 1, CalculateGeneResponse, sex_key=females_key)
response_females_x <- apply(expr_females_x, 1, CalculateGeneResponse, sex_key=females_key)
response_males_auto <- apply(expr_males_auto, 1, CalculateGeneResponse, sex_key=males_key)
response_males_x <- apply(expr_males_x, 1, CalculateGeneResponse, sex_key=males_key)

#Make combined dataframes
foldchange_df <- MakeCombinedDataFrame(response_females_auto, response_males_auto, response_females_x, response_males_x)
foldchange_var_df <- MakeCombinedDataFrame(response_females_auto_var, response_males_auto_var, response_females_x_var, response_males_x_var)
foldchange_sampled_df <- MakeCombinedDataFrame(response_females_auto_final, response_males_auto_var, response_females_x_final, response_males_x_var)

#Sort dataframes
foldchange_df <- foldchange_df[order(foldchange_df$gene),]
foldchange_sampled_df <- foldchange_sampled_df[order(foldchange_sampled_df$gene),]
foldchange_var_df <- foldchange_var_df[order(foldchange_var_df$gene),]

#Make CV and Fano dataframes
CV_sampled_df <- cbind.data.frame(gene=foldchange_df$gene, chr_type=foldchange_df$chr_type, 
                       females_CV=abs(sqrt(foldchange_sampled_df$females_stat)/foldchange_df$females_stat),
                       males_CV=abs(sqrt(foldchange_sampled_df$males_stat)/foldchange_df$males_stat))
CV_df <- cbind.data.frame(gene=foldchange_df$gene, chr_type=foldchange_df$chr_type, 
                                  females_CV=abs(sqrt(foldchange_var_df$females_stat)/foldchange_df$females_stat),
                                  males_CV=abs(sqrt(foldchange_var_df$males_stat)/foldchange_df$males_stat))
Fano_sampled_df <- cbind.data.frame(gene=foldchange_df$gene, chr_type=foldchange_df$chr_type, 
                                  females_Fano=abs((foldchange_sampled_df$females_stat)/foldchange_df$females_stat),
                                  males_Fano=abs((foldchange_sampled_df$males_stat)/foldchange_df$males_stat))

#Calculate regression lines for sampled plot
auto_lm <- lm(CV_sampled_df[which(CV_sampled_df$chr_type == "auto"),"males_CV"] ~ CV_sampled_df[which(CV_sampled_df$chr_type == "auto"),"females_CV"] - 1)
x_lm <- lm(CV_sampled_df[which(CV_sampled_df$chr_type == "x"),"males_CV"] ~ CV_sampled_df[which(CV_sampled_df$chr_type == "x"),"females_CV"] - 1)
auto_Fano_lm <- lm(Fano_sampled_df[which(Fano_sampled_df$chr_type == "auto"),"males_Fano"] ~ Fano_sampled_df[which(Fano_sampled_df$chr_type == "auto"),"females_Fano"] - 1)
x_Fano_lm <- lm(Fano_sampled_df[which(Fano_sampled_df$chr_type == "x"),"males_Fano"] ~ Fano_sampled_df[which(Fano_sampled_df$chr_type == "x"),"females_Fano"] - 1)


#Plot Coefficients of Variances and Fano Factors
out_dir <- "C:/Users/Janel/Documents/sex-specific/Visuals/variance_plots/"
tiff(paste(out_dir, "CV_Comparison_sampled", ".jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
print(suppressWarnings(ggplot(CV_sampled_df, aes(x = females_CV, y = males_CV)) + 
                         geom_point(aes(col=chr_type)) + #add points colored by significance and nobody cares about escape genes
                         scale_color_manual(values=c(alpha("#E95C20FF", 0.4), "#00ff00")) +
                         ylim(min(CV_df$males_CV),max(CV_df$males_CV)) + 
                         xlim(min(CV_df$females_CV),max(CV_df$females_CV)) + 
                         #theme_classic() + 
                         theme(legend.position="none", legend.title = element_blank()) +
                         #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                         geom_hline(yintercept = 0, linetype="solid", color="black") + 
                         geom_vline(xintercept = 0, linetype="solid", color="black") + 
                         geom_abline(slope = 1, intercept = 0, size=0.1, linetype="dashed") + 
                         geom_abline(slope = x_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color="#00ff00") + #auto regression line
                         geom_abline(slope = auto_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("#E95C20FF", 1)) + #x regression line
                         #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                         xlab("Female Coefficient of Variation") + ylab("Male Coefficient of Variation")))
dev.off()

tiff(paste(out_dir, "Fano_Comparison_sampled", ".jpg", sep = ""), width = 8, height = 8, units = 'in', res = 500)
print(suppressWarnings(ggplot(Fano_sampled_df, aes(x = females_Fano, y = males_Fano)) + 
                         geom_point(aes(col=chr_type)) + #add points colored by significance and nobody cares about escape genes
                         scale_color_manual(values=c(alpha("#E95C20FF", 0.4), "#00ff00")) +
                         ylim(min(Fano_sampled_df$males_Fano),max(Fano_sampled_df$males_Fano)) + 
                         xlim(min(Fano_sampled_df$females_Fano),max(Fano_sampled_df$females_Fano)) + 
                         #theme_classic() + 
                         theme(legend.position="none", legend.title = element_blank()) +
                         #geom_text_repel(aes(label=Name, color = factor(input_add$Expression)), size=3, hjust=0, vjust=0, show.legend = FALSE, force = 3) + 
                         geom_hline(yintercept = 0, linetype="solid", color="black") + 
                         geom_vline(xintercept = 0, linetype="solid", color="black") + 
                         geom_abline(slope = 1, intercept = 0, size=0.1, linetype="dashed") + 
                         geom_abline(slope = x_Fano_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color="#00ff00") + #auto regression line
                         geom_abline(slope = auto_Fano_lm$coefficients[1], intercept = 0, size=0.5, linetype="dashed", color=alpha("#E95C20FF", 1)) + #x regression line
                         #annotate("text", x=1.5, y=-3, label = paste("correlation = ", format(cor(input_add$females_FC, input_add$males_FC), digits=3), sep = ""), hjust = -0.05, vjust = 1) +
                         xlab("Female Fano Factor") + ylab("Male Fano Factor")))
dev.off()

#Run autosomal Wilcox Test on Fano Factors and CVs
Fano_autosomal_test <- wilcox.test(Fano_sampled_df[which(Fano_sampled_df$chr_type == "auto"),"females_Fano"], 
                                      Fano_sampled_df[which(Fano_sampled_df$chr_type == "auto"),"males_Fano"], paired = TRUE, 
                                      alternative = "greater")
CV_autosomal_test <- wilcox.test(CV_sampled_df[which(CV_sampled_df$chr_type == "auto"),"females_CV"], 
                                   CV_sampled_df[which(CV_sampled_df$chr_type == "auto"),"males_CV"], paired = TRUE, 
                                   alternative = "two.sided")

# 13) Execute variance tests ====

#Sort all data frames to align rows
#expr_females_x <- expr_females_x[order(rownames(expr_females_x)),]
#expr_males_x <- expr_males_x[order(rownames(expr_males_x)),]
#expr_females_auto <- expr_females_auto[order(rownames(expr_females_auto)),]
#expr_males_auto <- expr_males_auto[order(rownames(expr_males_auto)),]

#Create lists of results needed for loops
#auto_LPS <- vector()
#x_LPS <- vector()
#auto_VEH <- vector()
#x_VEH <- vector()

#Cycle through rows of x matrices first
#for (i in 1:nrow(expr_females_x)) {
  
#  #Create dfs of data needed for gene (by gender)
#  males_df <- as.data.frame(cbind(extract_ord=(colnames(expr_males_x)), gene_expr=as.numeric(t(expr_males_x[i,])), gender=(rep("males", ncol(expr_males_x)))))
#  females_df <- as.data.frame(cbind(extract_ord=(colnames(expr_females_x)), gene_expr=as.numeric(t(expr_females_x[i,])), gender=(rep("females", ncol(expr_females_x)))))
  
#  #split gender dfs into separate gender/treatment categories
#  LPS_df <- rbind(males_df[which(males_df$extract_ord %in% paste("X", males_key[,"LPS"], sep = "")),], females_df[which(females_df$extract_ord %in% paste("X", females_key[,"LPS"], sep = "")),])
#  VEH_df <- rbind(males_df[which(males_df$extract_ord %in% paste("X", males_key[,"null"], sep = "")),], females_df[which(females_df$extract_ord %in% paste("X", females_key[,"null"], sep = "")),])
  
#  #Make sure gene_expr columns are numeric
#  LPS_df$gene_expr <- as.numeric(as.character(LPS_df$gene_expr))
#  VEH_df$gene_expr <- as.numeric(as.character(VEH_df$gene_expr))
  
  #Run LPS test
#  LPS_test <- bf.test(gene_expr ~ gender, LPS_df)
#  VEH_test <- bf.test(gene_expr ~ gender, VEH_df)
#  
#  #Append to results vectors
#  x_LPS <- append(x_LPS, LPS_test$p.value)
#  x_VEH <- append(x_VEH, VEH_test$p.value)
#}

#Cycle through rows of auto matrices second
#for (i in 1:nrow(expr_females_auto)) {
#  
#  #Create dfs of data needed for gene (by gender)
#  males_df <- as.data.frame(cbind(extract_ord=(colnames(expr_males_auto)), gene_expr=as.numeric(t(expr_males_auto[i,])), gender=(rep("males", ncol(expr_males_auto)))))
#  females_df <- as.data.frame(cbind(extract_ord=(colnames(expr_females_auto)), gene_expr=as.numeric(t(expr_females_auto[i,])), gender=(rep("females", ncol(expr_females_auto)))))
#  
#  #split gender dfs into separate gender/treatment categories
#  LPS_df <- rbind(males_df[which(males_df$extract_ord %in% paste("X", males_key[,"LPS"], sep = "")),], females_df[which(females_df$extract_ord %in% paste("X", females_key[,"LPS"], sep = "")),])
#  VEH_df <- rbind(males_df[which(males_df$extract_ord %in% paste("X", males_key[,"null"], sep = "")),], females_df[which(females_df$extract_ord %in% paste("X", females_key[,"null"], sep = "")),])
#  
#  #Make sure gene_expr columns are numeric
#  LPS_df$gene_expr <- as.numeric(as.character(LPS_df$gene_expr))
#  VEH_df$gene_expr <- as.numeric(as.character(VEH_df$gene_expr))
#  
#  #Run LPS test
#  LPS_test <- bf.test(gene_expr ~ gender, LPS_df)
#  VEH_test <- bf.test(gene_expr ~ gender, VEH_df)
#  
#  #Append to results vectors
#  auto_LPS <- append(auto_LPS, LPS_test$p.value)
#  auto_VEH <- append(auto_VEH, VEH_test$p.value)
#}

# 14) Correct for Multiple Testing and analyze results ====

#First remove junk
#remove(i, filter_dir, females_df, males_df, LPS_df, VEH_df, VEH_test, LPS_test)
#remove(inp_dir)

#Adjust for multiple testing??????????????????????????????????????????
#auto_LPS_corrected <- p.adjust(auto_LPS, method = "fdr")
#auto_VEH_corrected <- p.adjust(auto_VEH, method = "fdr")
#x_LPS_corrected <- p.adjust(x_LPS, method = "fdr")
#x_VEH_corrected <- p.adjust(x_VEH, method = "fdr")

#Make count variables
#counts_auto_LPS <- c(length(auto_LPS_corrected[which(auto_LPS_corrected <= 0.05)]),length(auto_LPS_corrected)-length(auto_LPS_corrected[which(auto_LPS_corrected <= 0.05)]))
#counts_auto_VEH <- c(length(auto_VEH_corrected[which(auto_VEH_corrected <= 0.05)]),length(auto_VEH_corrected)-length(auto_VEH_corrected[which(auto_VEH_corrected <= 0.05)]))
#counts_x_LPS <- c(length(x_LPS_corrected[which(x_LPS_corrected <= 0.05)]),length(x_LPS_corrected)-length(x_LPS_corrected[which(x_LPS_corrected <= 0.05)]))
#counts_x_VEH <- c(length(x_VEH_corrected[which(x_VEH_corrected <= 0.05)]),length(x_VEH_corrected)-length(x_VEH_corrected[which(x_VEH_corrected <= 0.05)]))

#Make test tables
#LPS_table <- cbind.data.frame(counts_auto_LPS, counts_x_LPS)
#VEH_table <- cbind.data.frame(counts_auto_VEH, counts_x_VEH)
#auto_table <- cbind.data.frame(counts_auto_LPS, counts_auto_VEH) 
#x_table <- cbind.data.frame(counts_x_LPS, counts_x_VEH) 

#Execute fisher tests
#p.values <- c("LPS"=fisher.test(LPS_table)$p.value, "VEH"=fisher.test(VEH_table)$p.value, "auto"=fisher.test(auto_table)$p.value, "x"=fisher.test(x_table)$p.value)
#p.values <- p.adjust(p.values, method = "fdr")


