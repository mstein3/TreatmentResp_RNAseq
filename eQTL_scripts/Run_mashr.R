##### 1. Set-up terminal on Gardner ####
# Start: 1.23.20


#Call needed libraries (They are already installed on Gardner)
library(assertthat)
library(mvtnorm)
library(rmeta)
library(profmem)
library(REBayes)
library(mashr)
#Set directories
setwd("/scratch/mconery/sex-specific/mashr")
inp_dir <- "/scratch/mconery/sex-specific/eQTL_analysis"
#Testing inp_dir
#inp_dir <- "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/eQTL_analysis/raw_data"

##### 2. Organize the files to get the SNPs with the max z score (Run Offline) ==========

#Read in file
mashr_raw <- read.delim(file = paste(inp_dir, "combined.auto.gemma.assoc.tsv", sep = "/"), 
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#Set total tests and clean up colnames of df
total_tests <- as.numeric(mashr_raw[1,1])
colnames(mashr_raw) <- mashr_raw[2,]
mashr_raw <- mashr_raw[3:nrow(mashr_raw),]

#Convert Min_p_score to numeric
mashr_raw$Min_p_score <- as.numeric(mashr_raw$Min_p_score)

#sort and then subset out strongest results for each gene
mashr_raw <- mashr_raw[order(mashr_raw$Min_p_score),]
mashr_strong <- mashr_raw[which(duplicated(mashr_raw$gene) == FALSE),]

#Resort dataframes
mashr_raw <- mashr_raw[order(mashr_raw$gene),]
mashr_strong <- mashr_strong[order(mashr_strong$gene),]

#Recharacterize columns as numerics
mashr_strong$ML_beta <- as.numeric(mashr_strong$ML_beta)
mashr_strong$MV_beta <- as.numeric(mashr_strong$MV_beta)
mashr_strong$FL_beta <- as.numeric(mashr_strong$FL_beta)
mashr_strong$FV_beta <- as.numeric(mashr_strong$FV_beta)
mashr_raw$ML_beta <- as.numeric(mashr_raw$ML_beta)
mashr_raw$MV_beta <- as.numeric(mashr_raw$MV_beta)
mashr_raw$FL_beta <- as.numeric(mashr_raw$FL_beta)
mashr_raw$FV_beta <- as.numeric(mashr_raw$FV_beta)
mashr_strong$ML_se <- as.numeric(mashr_strong$ML_se)
mashr_strong$MV_se <- as.numeric(mashr_strong$MV_se)
mashr_strong$FL_se <- as.numeric(mashr_strong$FL_se)
mashr_strong$FV_se <- as.numeric(mashr_strong$FV_se)
mashr_raw$ML_se <- as.numeric(mashr_raw$ML_se)
mashr_raw$MV_se <- as.numeric(mashr_raw$MV_se)
mashr_raw$FL_se <- as.numeric(mashr_raw$FL_se)
mashr_raw$FV_se <- as.numeric(mashr_raw$FV_se)

# make dataframes of se and beta for "strong" data
strong_se<-cbind.data.frame(MV_se=mashr_strong$MV_se, ML_se=mashr_strong$ML_se, FV_se=mashr_strong$FV_se, FL_se=mashr_strong$FL_se)
strong_beta<-cbind.data.frame(MV_beta=mashr_strong$MV_beta, ML_beta=mashr_strong$ML_beta, FV_beta=mashr_strong$FV_beta, FL_beta=mashr_strong$FL_beta)
rownames(strong_se)<-paste(mashr_strong$gene, mashr_strong$rs, sep = ":")
rownames(strong_beta)<-paste(mashr_strong$gene, mashr_strong$rs, sep = ":")
colnames(strong_se)<-c("MaleNull", "MaleLPS", "FemaleNull", "FemaleLPS")
colnames(strong_beta)<-c("MaleNull", "MaleLPS", "FemaleNull", "FemaleLPS")

# make dataframes of se and beta for "random" data
random_se<-cbind.data.frame(MV_se=mashr_raw$MV_se, ML_se=mashr_raw$ML_se, FV_se=mashr_raw$FV_se, FL_se=mashr_raw$FL_se)
random_beta<-cbind.data.frame(MV_beta=mashr_raw$MV_beta, ML_beta=mashr_raw$ML_beta, FV_beta=mashr_raw$FV_beta, FL_beta=mashr_raw$FL_beta)
rownames(random_se)<-paste(mashr_raw$gene, mashr_raw$rs, sep = ":")
rownames(random_beta)<-paste(mashr_raw$gene, mashr_raw$rs, sep = ":")
colnames(random_se)<-c("MaleNull", "MaleLPS", "FemaleNull", "FemaleLPS")
colnames(random_beta)<-c("MaleNull", "MaleLPS", "FemaleNull", "FemaleLPS")


####### 3. Create the beta and SE files: (Run On Gardner) =========

#Remove existing junk
suppressWarnings(remove(mashr_raw, mashr_strong))

#Create list
strong_list <- list(beta_strong = data.matrix(strong_beta), se_strong = data.matrix(strong_se))
random_list <- list(beta_random = data.matrix(random_beta), se_random = data.matrix(random_se))

####### 4. Start running some of mash (Run On Gardner): =============
##### 4.1 Set up basic data-driven covariances ####
# from https://stephenslab.github.io/mashr/articles/eQTL_outline.html

#Read all data into mash format to reference in correlation structure creation
#This has been changed from previous run
data.temp <- mash_set_data(random_list$beta_random, random_list$se_random)

#Create correlation structure
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

#Set up main data objects
data.random<-mash_set_data(random_list$beta_random, random_list$se_random, V=Vhat)
data.strong<-mash_set_data(strong_list$beta_strong, strong_list$se_strong, V=Vhat)

#Set-up data driven covariances
U.pca = cov_pca(data.strong,4)
U.ed = cov_ed(data.strong, U.pca)

##### 4.2 Fit mash model (estimate mixture proportions) ######
#Now we fit mash to the random tests using both data-driven and canonical covariances. 
#(Remember the Crucial Rule! We have to fit using a random set of tests, and not a dataset that is enriched for strong tests.) 
#The outputlevel=1 option means that it will not compute posterior summaries for these tests (which saves time).

#Changed again to random tests
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

#Here are the console outputs
#- Computing 848116 x 393 likelihood matrix.
#- Likelihood calculations took 526.86 seconds.
#- Fitting model with 393 mixture components.


##### 4.3 Compute Posterior Summaries ######
#Now we can compute posterior summaries etc for any subset of tests using the above mash fit. 
#Here we do this for the strong tests. We do this using the same mash function as above, 
#but we specify to use the fit from the previous run of mash by specifying
#g=get_fitted_g(m), fixg=TRUE. (In mash the parameter g is used to denote the mixture model which we learned above.)

#Run for strong data only
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
#Console outputs


head(get_lfsr(m2))
#Console outputs


#Set wd to results folder
out_dir <- "/scratch/mconery/sex-specific/mashr/"
saveRDS(m2, paste(out_dir, "Mash_Results_allgenes_strong_020320.rds", sep = ""))

#Get significant results and p-values
lfsr_results<-get_lfsr(m2)
#Output p-value results
write.table(lfsr_results, paste(out_dir, "AllRes_lfsr_combined_031120.txt", sep = ""), sep="\t", quote=FALSE)

##Get significant results
head(get_significant_results(m2))
print(length(get_significant_results(m2)))
all_sig<-get_significant_results(m2)

#Can also get the significant results in a subset of conditions. 
print(length(get_significant_results(m2, conditions="MaleNull")))
print(length(get_significant_results(m2, conditions="MaleLPS")))
print(length(get_significant_results(m2, conditions="FemaleNull")))
print(length(get_significant_results(m2, conditions="FemaleLPS")))
