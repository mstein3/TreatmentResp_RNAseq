#Testing data
#GRM_dir <- "C:/Users/Janel/Documents/sex-specific/prep_GRMs/"
#process_dir <- "C:/Users/Janel/Documents/sex-specific/processed_data/"
#filter_dir <- "C:/Users/Janel/Documents/sex-specific/filters/"


#Internal function that takes a restricted grm and makes a final grm using the pheno ids as inputs
#Needs a pheno_file pre-sorted by Extract_ords 
make_grm <- function(grm_int, pheno_raw){
  
  #Create blank vector to receive values
  final_grm_vector <- vector()
  
  #Cycle through extract_ord ids in pheno file and make rows for all ids below the row in the file
  for (i in 1:nrow(pheno_raw)) {
    #Check that FINDIV corresponding to id is in the grm_int (some findivs are not)
    if (pheno_raw[i,"FINDIV"] %in% grm_int[,1]) {
      #Cycle through rows in pheno df from i to bottom of matrix (need self value)
      #Originally this was i:nrow(pheno_raw)
      for (j in i:nrow(pheno_raw)) {
        #Verify that FINDIV of j row is also in grm_int
        if (pheno_raw[j,"FINDIV"] %in% grm_int[,2]) {
          #Check which FINDIV is smaller to know order of columns in grm_int
          #Put in 10000 as a place holder
          if (pheno_raw[i,"FINDIV"] <= pheno_raw[j,"FINDIV"]) {
            final_grm_vector <- append(final_grm_vector, paste(pheno_raw[i,"Extract_Ord"], pheno_raw[j,"Extract_Ord"], 
                                                               grm_int[which(grm_int$V1==pheno_raw[i,"FINDIV"] & 
                                                                               grm_int$V2==pheno_raw[j,"FINDIV"]),"V3"], sep = ","))
          } else {
            final_grm_vector <- append(final_grm_vector, paste(pheno_raw[i,"Extract_Ord"], pheno_raw[j,"Extract_Ord"], 
                                                               grm_int[which(grm_int$V1==pheno_raw[j,"FINDIV"] & 
                                                                               grm_int$V2==pheno_raw[i,"FINDIV"]),"V3"], sep = ","))
          }
        }
      }
    }
  }
  
  #Separate vector into data frame
  final_grm <- data.frame(matrix(unlist(strsplit(final_grm_vector, split = ",")), nrow = length(final_grm_vector), ncol = 3, byrow = TRUE))
  #Sort data frame by extract_ords before returning
  final_grm$X1 <- as.numeric(as.character(final_grm$X1))
  final_grm$X2 <- as.numeric(as.character(final_grm$X2))
  final_grm  <- final_grm[order(final_grm$X1, final_grm$X2),]
  #return data frame
  return(final_grm)
}

#Outer function that splits grm by gender, calls internal function, and then writes to file
clean_grm <- function(grm_raw, females_pheno, males_pheno, write_loc, filter_write_loc){
  
  #Combine data from pheno files
  combined_pheno <- rbind(females_pheno[,2:7], males_pheno[,2:7])
  #Need to sort combined_pheno by FINDIVs to ensure match happens correctly
  combined_pheno <- combined_pheno[order(combined_pheno$FINDIV),]
  #Create key of ids
  individuals <- combined_pheno[,"FINDIV"]
  individuals <- subset(individuals,duplicated(individuals)!=TRUE)
  combined_pheno$treat <- as.factor(combined_pheno$treat)
  id_key <- cbind.data.frame(individuals, combined_pheno[which(combined_pheno$treat=="LPS" & combined_pheno$FINDIV%in%individuals),"Extract_Ord"], 
                  combined_pheno[which(combined_pheno$treat=="null" & combined_pheno$FINDIV%in%individuals),"Extract_Ord"],
                  combined_pheno[which(combined_pheno$treat=="null" & combined_pheno$FINDIV%in%individuals),"sex"])
  colnames(id_key) <-  c("FINDIV", "LPS", "VEH", "SEX")
  remove(individuals)

  
  #remove junk rows and split grms by gender
  grm_int <- grm_raw[which(grm_raw$V1 %in% id_key[,"FINDIV"] & grm_raw$V2 %in% id_key[,"FINDIV"]),]
  remove(grm_raw)
  
  #Call make_grm function
  final_grm <- make_grm(grm_int, combined_pheno)
  #Remove grm_int
  remove(grm_int)
  
  #Write to files (tab delimited)
  write.table(final_grm, file = paste(write_loc, "overall.gemma.grm", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  #Create id file 
  combined_id <- cbind.data.frame(fid=rep("HUTTERITES", nrow(combined_pheno)), iid=combined_pheno[,"Extract_Ord"])
  combined_id <- combined_id[order(combined_id$iid),]
  write.table(combined_id, file = paste(filter_write_loc, "overall.iids.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}


# Define input directories
GRM_dir <- "/scratch/mconery/sex-specific/GRMs/"
process_dir <- "/scratch/mconery/sex-specific/processed_data/"
filter_dir <- "/scratch/mconery/sex-specific/filters/"

#Read in pheno files
females_pheno <- read.table(paste(process_dir, "females.Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
males_pheno <- read.table(paste(process_dir, "males.Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Sort pheno files
attach(females_pheno)
females_pheno <- females_pheno[order(Extract_Ord),]
detach(females_pheno)
attach(males_pheno)
males_pheno <- males_pheno[order(Extract_Ord),]
detach(males_pheno)

#Read in filter files and filter out FINDIV not in ids lists
females_filter <- read.table(paste(filter_dir, "females.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
males_filter <- read.table(paste(filter_dir, "males.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
females_pheno <- females_pheno[which(females_pheno$FINDIV %in% females_filter$V2),]
males_pheno <- males_pheno[which(males_pheno$FINDIV %in% males_filter$V2),]
remove(females_filter, males_filter)


#read in raw auto-grm and clean it
auto_grm_raw <- read.table(paste(GRM_dir,"additive.coef.2015-12-15", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
clean_grm(auto_grm_raw, females_pheno, males_pheno, GRM_dir, filter_dir)
remove(auto_grm_raw)

