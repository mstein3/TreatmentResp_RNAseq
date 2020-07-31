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
  #return data frame
  return(final_grm)
}

#Outer function that splits grm by gender, calls internal function, and then writes to file
clean_grm <- function(grm_raw, females_pheno, males_pheno, write_loc, filter_write_loc, grm_type, treatment){
  
  #Combine data from pheno files
  combined_pheno <- rbind(females_pheno[,2:7], males_pheno[,2:7])
  #Need to sort combined_pheno by FINDIVs first (will sort back later)
  combined_pheno <- combined_pheno[order(combined_pheno$FINDIV),]
  #Create key of ids
  individuals <- combined_pheno[,"FINDIV"]
  individuals <- subset(individuals,duplicated(individuals)!=TRUE)
  combined_pheno$treat <- as.factor(combined_pheno$treat)
  #Check for treatment type to know which column extract ord to pull
  #Also filter down females and males_pheno files for just the treatment type
  if (treatment == "treatment") {
    id_key <- cbind.data.frame(individuals, combined_pheno[which(combined_pheno$treat=="LPS" & combined_pheno$FINDIV%in%individuals),"Extract_Ord"], 
                    combined_pheno[which(combined_pheno$treat=="LPS" & combined_pheno$FINDIV%in%individuals),"sex"])
    females_pheno <- females_pheno[which(females_pheno$treat == "LPS"),]
    males_pheno <- males_pheno[which(males_pheno$treat == "LPS"),]
  } else if (treatment == "no.treatment") {
    id_key <- cbind.data.frame(individuals, combined_pheno[which(combined_pheno$treat=="null" & combined_pheno$FINDIV%in%individuals),"Extract_Ord"], 
                    combined_pheno[which(combined_pheno$treat=="null" & combined_pheno$FINDIV%in%individuals),"sex"])
    females_pheno <- females_pheno[which(females_pheno$treat == "null"),]
    males_pheno <- males_pheno[which(males_pheno$treat == "null"),]
  }
  colnames(id_key) <-  c("FINDIV", "Extract_Ord", "SEX")
  remove(individuals, combined_pheno)
  #Resort id_key by extract_ords
  id_key <- id_key[order(id_key$Extract_Ord),]
  
  #remove junk rows and split grms by gender
  grm_females_int <- grm_raw[which(grm_raw$V1 %in% id_key[which(id_key[,"SEX"]==1),"FINDIV"] & grm_raw$V2 %in% id_key[which(id_key[,"SEX"]==1),"FINDIV"]),]
  grm_males_int <- grm_raw[which(grm_raw$V1 %in% id_key[which(id_key[,"SEX"]==0),"FINDIV"] & grm_raw$V2 %in% id_key[which(id_key[,"SEX"]==0),"FINDIV"]),]
  remove(grm_raw)
  
  #Call make_grm function
  females_grm <- make_grm(grm_females_int, females_pheno)
  males_grm <- make_grm(grm_males_int, males_pheno)
  #Remove grm_int
  remove(grm_males_int, grm_females_int)
  
  #Write to files (tab delimited)
  write.table(females_grm, file = paste(write_loc, "females.", grm_type, ".", treatment, ".new.gemma.grm", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(males_grm, file = paste(write_loc, "males.", grm_type, ".", treatment, ".new.gemma.grm", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  #Create id files comment out
  females_id <- females_pheno[which(females_pheno$Extract_Ord %in% females_grm[,1]),]
  males_id <- males_pheno[which(males_pheno$Extract_Ord %in% males_grm[,1]),]
  females_id <- cbind(females_id[,"FID"], females_id[,"Extract_Ord"])
  males_id <- cbind(males_id[,"FID"], males_id[,"Extract_Ord"])
  write.table(females_id, file = paste(filter_write_loc, "females.", grm_type, ".", treatment, ".iids.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(males_id, file = paste(filter_write_loc, "males.", grm_type, ".", treatment , ".iids.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}


# Define input directories
GRM_dir <- "/scratch/mconery/sex-specific/GRMs/"
process_dir <- "/scratch/mconery/sex-specific/processed_data/"
filter_dir <- "/scratch/mconery/sex-specific/filters/"

#Read in pheno files
females_pheno_outside <- read.table(paste(process_dir, "females.Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
males_pheno_outside <- read.table(paste(process_dir, "males.Pheno.Rready.nullLPS.7M.031520.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Sort pheno files
attach(females_pheno_outside)
females_pheno_outside <- females_pheno_outside[order(Extract_Ord),]
detach(females_pheno_outside)
attach(males_pheno_outside)
males_pheno_outside <- males_pheno_outside[order(Extract_Ord),]
detach(males_pheno_outside)

#Read in filter files and filter out FINDIV not in ids lists
females_filter <- read.table(paste(filter_dir, "females.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
males_filter <- read.table(paste(filter_dir, "males.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
females_pheno_outside <- females_pheno_outside[which(females_pheno_outside$FINDIV %in% females_filter$V2),]
males_pheno_outside <- males_pheno_outside[which(males_pheno_outside$FINDIV %in% males_filter$V2),]
remove(females_filter, males_filter)

#Define treatments and cycle through both treatment options
treatments <- c("treatment", "no.treatment")
for (treatment in treatments) {
  #Read in X GRM
  X_grm_raw <- read.table(paste(GRM_dir,"additive.X.coef.3555", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  #Read in filters
  x_grm_ids_females <- read.table(paste(filter_dir, "females.X.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  x_grm_ids_males <- read.table(paste(filter_dir, "males.X.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  #Filter pheno_files
  females_pheno_x <- females_pheno_outside[which(females_pheno_outside$Extract_Ord %in% x_grm_ids_females[,2]),]
  males_pheno_x <- males_pheno_outside[which(males_pheno_outside$Extract_Ord %in% x_grm_ids_males[,2]),]
  #Clean GRMs and write to file
  clean_grm(X_grm_raw, females_pheno_x, males_pheno_x, GRM_dir, filter_dir, "X", treatment)
  #Clean-up
  remove(X_grm_raw, x_grm_ids_females, x_grm_ids_males)
  
  #Read in auto GRM
  auto_grm_raw <- read.table(paste(GRM_dir,"additive.coef.2015-12-15", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  #Read in filters
  auto_grm_ids_females <- read.table(paste(filter_dir, "females.auto.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  auto_grm_ids_males <- read.table(paste(filter_dir, "males.auto.iids.txt", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  #Filter pheno_files
  females_pheno_auto <- females_pheno_outside[which(females_pheno_outside$Extract_Ord %in% auto_grm_ids_females[,2]),]
  males_pheno_auto <- males_pheno_outside[which(males_pheno_outside$Extract_Ord %in% auto_grm_ids_males[,2]),]
  #Clean GRMs and write to file
  clean_grm(auto_grm_raw, females_pheno_auto, males_pheno_auto, GRM_dir, filter_dir, "auto", treatment)
  #Clean-up
  remove(auto_grm_raw, auto_grm_ids_females, auto_grm_ids_males)
}


