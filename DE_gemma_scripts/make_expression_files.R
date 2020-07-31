# 0) Install and load libraries ====

#call libraries
#library(plyr)
#library(dplyr)
#library(tidyr)
#library(stringr)
#library(ggplot2)
#library(magrittr)
#library(scales)
#library(DescTools)


# 1) Read in files and match them by ids ====

#Cycle through genders and grm_types
genders = c("males", "females")
grm_types = c("auto", "X")

for (gender in genders) {
  for (grm_type in grm_types) {
    #Set input directory
    inp_dir <- "/scratch/mconery/sex-specific/processed_data/"
    filter_dir <- "/scratch/mconery/sex-specific/filters/"
    
    #Read in pheno file
    pheno_raw <- read.table(paste(inp_dir, gender, ".Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    #Read in filter files and filter out FINDIV not in ids list
    ids_filter <- read.table(paste(filter_dir, gender, ".", grm_type, ".iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    pheno_raw <- pheno_raw[which(pheno_raw$Extract_Ord %in% ids_filter$V2),]
    remove(ids_filter)
    #Sort pheno_file
    pheno_raw <- pheno_raw[order(pheno_raw$Extract_Ord),]
    
    #Check for GRM type and only read in auto file if necessary
    #X-grm set needed for both autosomal and x-linked genes
    if (grm_type == "auto") {
      #Read in expr files
      x_expr_raw <- read.table(paste(inp_dir, gender, ".logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = ""), sep = "\t", header = TRUE)
      auto_expr_raw <- read.table(paste(inp_dir, "BatchCorr_TMMvoom_NullLPS_GX_224_20190223.txt", sep = ""), sep = "\t", header = TRUE)
      #filter expr files
      #Use sprintf to make all numbers four digits, i.e. add leading zeroes if necessary
      x_expr_raw <- x_expr_raw[,sprintf("X%04d", pheno_raw$Extract_Ord)]
      auto_expr_raw <- auto_expr_raw[,paste("X", pheno_raw$Extract_Ord, sep="")]
      #Overwrite column names of x_expr_raw to match the format of auto_expr_raw
      colnames(x_expr_raw) <- colnames(auto_expr_raw)
      #Rbind x_expr_raw to auto_expr_raw to make expr_raw
      expr_raw <- rbind.data.frame(auto_expr_raw, x_expr_raw)
    } else { #If we are looking at x-linked data than only read in x-linked expression matrix
      expr_raw <- read.table(paste(inp_dir, gender, ".logCPM.NullLPS.RunsCombined.GeneNames.NoPoolExtract.031520.txt", sep = ""), sep = "\t", header = TRUE)
      expr_raw <- expr_raw[,sprintf("X%04d", pheno_raw$Extract_Ord)]
      colnames(expr_raw) <- paste("X", pheno_raw$Extract_Ord, sep = "")
    }

    
    #Get treatment/non-treatment iids 
    treatment_ids <- subset(pheno_raw, treat == "LPS")["Extract_Ord"]
    null_ids <- subset(pheno_raw, treat == "null")["Extract_Ord"]
    
    #write expression files
    for (i in 1:nrow(expr_raw)) {
      #create dfs
      pheno <- as.data.frame(cbind(pheno_raw[,"FID"], pheno_raw[,"Extract_Ord"], t(expr_raw[i,])))
      colnames(pheno) <- c("FID", "IID", rownames(expr_raw)[i])
      pheno_treatment <- subset(pheno, IID %in% treatment_ids[,1])
      pheno_no_treatment <- subset(pheno, IID %in% null_ids[,1])
      
      #write files
      write.table(pheno, file = paste("/scratch/mconery/sex-specific/expression_data/", gender, ".", rownames(expr_raw)[i], ".", grm_type, ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      write.table(pheno_treatment, file = paste("/scratch/mconery/sex-specific/expression_data/", gender, ".", rownames(expr_raw)[i], ".", grm_type, ".treatment.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      write.table(pheno_no_treatment, file = paste("/scratch/mconery/sex-specific/expression_data/", gender, ".", rownames(expr_raw)[i], ".", grm_type, ".no.treatment.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    } 
  }
}
