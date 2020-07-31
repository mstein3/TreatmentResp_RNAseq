#### Processing RNAseq gene counts after mapping with STAR
# Start: 8/8/17

####### 1. Count the number of mapped reads =============
setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/mapped_reads_outputs_hg19/")
file_list <- list.files()
map_read_counts<-data.frame(matrix(NA, nrow = length(file_list), ncol = 5))

for(i in 1:length(file_list)) {
    dataset <- read.table(file_list[i], sep="\t")
    total.reads<-sum(dataset[,2])
    unmap.reads<-sum(dataset[1:4,2])
    map.reads<- (total.reads-unmap.reads)
    percent.map<- (map.reads/total.reads)
    unmap.reads.strict<-dataset[1,2]
    map.reads.strict<- (total.reads-unmap.reads.strict)
    percent.map.strict<-(map.reads.strict/total.reads)
    map_read_counts[i,1]<-paste(file_list[i])
    map_read_counts[i,2]<-map.reads
    map_read_counts[i,3]<-percent.map
    map_read_counts[i,4]<-map.reads.strict
    map_read_counts[i,5]<-percent.map.strict
  
}

map.table.names<-c("name", "NumMapReads_allAmbig", "PercentMapReads_allAmbig", "NumMapReads_UnmapOnly", "PercentMapReads_UnmapOnly")
colnames(map_read_counts)<-map.table.names

setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/")
write.table(map_read_counts, "Mapped_read_counts_8817.txt", sep="\t", quote=FALSE, row.names=FALSE)
map_read_counts<-read.table("Mapped_read_counts_8817.txt", sep="\t", header=TRUE)


####### 2. Read in the gene count files for processing. ========
setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/mapped_reads_outputs_hg19/")

file_list_test<-file_list[1:5]

for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, sep="\t")
    colnames(temp_dataset)<-c("V1", paste(file), "V3", "V4")
    dataset<-cbind.data.frame(dataset, temp_dataset[,2])
    rm(temp_dataset)
  }
  
}

dim(dataset)
# [1] 57824  1372

### Remove columns 2, 3 and 4 (from the first file read in )
dataset2<-dataset[,-2:-4]

# Read in file list with unique extract ord names, all with 4 digits. 
# Name columns by cleaning up some of the file list: 
file_list4<-read.table("file_list4.txt", sep="\t", header=FALSE)
colnames(dataset2)<-file_list4$V1

# Write this file 
setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/")
write.table(dataset2, "GeneCounts_SamplesNotMerged_9117.txt", sep="\t", row.names=F)
dataset2<-read.table("GeneCounts_SamplesNotMerged_9117.txt", sep="\t", header=TRUE, row.names=1)

# Remove the first 4 rows from the dataset, which do not correspond to genes:
dataset4<-dataset2[-1:-4,]

# Total up number of reads for each seq run: 
each.run.map<-colSums(dataset4)

###### 3.  Add together the technical replicates =============
#unique(sapply(strsplit(names(dataset4), "_", fixed=TRUE), "[", 2) )
#[1] "C7KFNACXX" "c7l0hacxx" "C7ERJACXX" "C7K36ACXX" "C7KRWACXX" "C7KV3ACXX"
#[7] "C7WPYACXX" "C7KVUACXX" "c8294acxx" "GTGAAA"    "C6TLUACXX" "C8TYAACXX"
#[13] "CGTACG"    "GTTTCG"    "ATTCCT"    "TTAGGC"    "GATCAG"    "CGATGT"   
#[19] "ACAGTG"    "CCGTCC"    "ACTGAT"    "AGTCAA"    "ATCACG"    "AGTTCC"   
#[25] "CAGATC"    "GTGGCC"    "ATGTCA"    "C8YATACXX" "GAGTGG"    "ACTTGA"   
#[31] "CTTGTA"    "GCCAAT"    "TAGCTT"    "GTCCGC"    "GGCTAC"    

#unique(sapply(strsplit(names(dataset4), "_", fixed=TRUE), "[", 1) )
#[1] "1001" "1002" "1003" "1004" "1005" "1007" "1008" "1009" "1010" "1011"
#[11] "1012" "1013" "1014" "1015" "1016" "1017" "1018" "1019" "102"  "1020"
# through the end of the extract ORr numbers 

colnames_for_org <- unique(sapply(strsplit(names(dataset4), "_", fixed=TRUE), "[", 1) )
# Add together columns with matching extract ords, using the uniform 4-digit codes (ex: 0007, 0087, etc.)
res<-sapply(colnames_for_org, function(x) rowSums(dataset4 [, grep(paste(x,"_",sep=""), names(dataset4)), drop=FALSE] )  )
dim(res)

#### But, are samples like 7, 9 have way too many reads? 
uniq.map.reads<-colSums(res)
summary(uniq.map.reads)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#14070   8564000  10530000  16600000  15270000 241900000 
# Yaaaas this is correct, finally. 

# Write out combined gene count table (this table still has gene names as ENSG though)
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R")
write.table(res, "RawGeneCounts_RunsCombined_ENSG_names_noQC_9117.txt", sep="\t", quote=FALSE)
res<-read.table("RawGeneCounts_RunsCombined_ENSG_names_noQC_9117.txt", sep="\t", header=TRUE)

#### Check: how many reads that map uniquely to genes: 
uniq.map.reads<-colSums(res)
reads.10M<-which(uniq.map.reads>10000000)
length(reads.10M)
# [1] 232 
reads.8M<-which(uniq.map.reads>8000000)
length(reads.8M)
# [1] 338
reads.7M<-which(uniq.map.reads>7000000)
length(reads.7M)
# [1] 382  
reads.less.3M<-which(uniq.map.reads<3000000)
length(reads.less.3M)
# [1] 6 (don't bother with these samples again)

####### 4. Remove X, Y, M genes ==================
setwd("/Users/kmagnaye/Documents/X_immunity_project")
ensg.all<-read.table("ENSG_GeneNames_Chr_All.txt", sep=" ", header=FALSE)

# Subset out protein coding genes:
protcode<-ensg.all$V5=="protein_coding"
ensg.pc<-ensg.all[protcode,]
dim(ensg.pc)
#[1] 20345     6

### 4.1 Keep only the X chr: 
x.chr<-ensg.pc$V1=="chrX"
ensg.x<-ensg.pc[x.chr,]
dim(ensg.x)
#[1] 830   6
length(unique(ensg.x$V6))
# [1] 829 (1 gene duplicated; IDS)

# Read in Gene Counts from raw data file and subset out protein-coding X chromosome genes from ensg file: 
res <- read.table("RawGeneCounts_RunsCombined_ENSG_names_noQC_9117.txt", sep = "\t", header = TRUE)
colnames(res) <- substr(colnames(res), 2, length(colnames(res)))
ensg.intersect<-intersect(rownames(res), ensg.x[,4])
length(ensg.intersect)
# [1] 830
ensg.match<-match(ensg.intersect, ensg.x[,4])
ensg.x2<-ensg.x[ensg.match,]
length(unique(ensg.x2$V6))
# [1] 829

# Subset: Subset out protein coding + autosomal genes + genes from res from res file: 
res.match<-match(ensg.intersect, rownames(res))
res2<-res[res.match,]
dim(res2)
#[1] 830 403

######### 5. Remove samples with <7M mapped reads: =========
pheno.start<-read.table("Pheno_Rready_403comb_8117.txt", header=TRUE, sep="\t")
enough.reads<-pheno.start$Total_Map>7000000
table(enough.reads)
# FALSE  TRUE 
# 24   379 
res3<-res2[,enough.reads]
dim(res3)
# [1] 830   379
pheno_7M<-pheno.start[enough.reads,]

######## 5.5 Subset by sex

females <- pheno_7M$sex==1
pheno_7M_female <- pheno_7M[females,"Extract_Ord"]
pheno_7M_male <- pheno_7M[!females,"Extract_Ord"]

res3_female <- res3[,colnames(res3) %in% pheno_7M_female]
res3_male <- res3[,colnames(res3) %in% pheno_7M_male]

######## 6. Pre-processing from the RNAseq analysis is as easy as 1-2-3 with limma, Glimma, and edgeR (subsetting lowly expressed genes)================
library(limma)
library(edgeR)

# Convert to CPM: 
cpm_female <- cpm(res3_female)
cpm_male <- cpm(res3_male)

lcpm_female <- cpm(res3_female, log=TRUE)
lcpm_male <- cpm(res3_male, log=TRUE)

# How many genes have a count of 0 in all samples? 
table(rowSums(res3_female==0)==199)
#FALSE  TRUE 
#666   164 

table(rowSums(res3_male==0)==121)
#FALSE  TRUE 
#660   170 

# How many genes have a count of 0 in both males and females?

female_0 <- rownames(res3_female)[rowSums(res3_female==0)==199]
male_0 <- rownames(res3_male)[rowSums(res3_male==0)==121]

both_0 <- intersect(female_0, male_0)
# intersect(female_0, male_0) == 145 genes with count 0 in both sexes
# female_0[!female_0 %in% both_0] == 19 genes with count 0 only in females
# male_0[!male_0 %in% both_0] == 25 genes with count 0 only in males

pdf("Log2CPM_female_dist_two_cutoffs_1117.pdf")
plot(density(lcpm_female[,12]), lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="A. Raw data in Females", xlab="Log-cpm")
abline(v=0, lty=3)
dev.off()

pdf("Log2CPM_male_dist_two_cutoffs_1117.pdf")
plot(density(lcpm_male[,1]), lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="A. Raw data in Males", xlab="Log-cpm")
abline(v=0, lty=3)
dev.off()

# With a typical cutoff. 
keep.exprs_female <- rowSums(cpm_female>1)>=10
cpm_female2 <- cpm_female[keep.exprs_female,]
dim(cpm_female2)
#[1] 548 199

keep.exprs_male <- rowSums(cpm_male>1)>=10
cpm_male2 <- cpm_male[keep.exprs_male,]
dim(cpm_male2)
#[1] 526 121

# Pull out from gene counts as well: 
cutoff.intersect_female<-intersect(rownames(cpm_female2), rownames(res3_female))
cutoff.match_female<-match(cutoff.intersect_female, rownames(res3_female))
res4_female<-res3_female[cutoff.match_female,]

cutoff.intersect_male<-intersect(rownames(cpm_male2), rownames(res3_male))
cutoff.match_male<-match(cutoff.intersect_male, rownames(res3_male))
res4_male<-res3_male[cutoff.match_male,]

# There are 525 genes in females and 501 in males
setwd("/Users/kmagnaye/Documents/X_immunity_project/Processing_in_R")
write.table(res4_female, "GeneCounts_Female_ENSG_noQC_12kgenes_110917.txt", sep="\t", quote=FALSE)
write.table(res4_male, "GeneCounts_Male_ENSG_noQC_12kgenes_110917.txt", sep="\t", quote=FALSE)
write.table(cpm_female2, "CPM_Female_ENSG_noQC_12kgenes_110917.txt", sep="\t", quote=FALSE)
write.table(cpm_male2, "CPM_Male_ENSG_noQC_12kgenes_110917.txt", sep="\t", quote=FALSE)

##### 7.  Try to convert to gene names============
ensg.intersect.12k_female<-intersect(rownames(res4_female), ensg.x[,4])
ensg.intersect.12k_male<-intersect(rownames(res4_male), ensg.x[,4])
#res.match.12k<-match(ensg.intersect.12k, rownames(res4))
ensg.match.12k_female<-match(ensg.intersect.12k_female, ensg.x[,4])
ensg.match.12k_male<-match(ensg.intersect.12k_male, ensg.x[,4])
ensg.12k_female<-ensg.x[ensg.match.12k_female,]
ensg.12k_male<-ensg.x[ensg.match.12k_male,]

length(unique(ensg.12k_female[,6]))
# [1] 547
length(unique(ensg.12k_male[,6]))
# [1] 525
# So, there's one dup in the females and males

# Re-order rownames in both data sets (res4, ensg12k)
ensg.12kb_female<-ensg.12k_female[match(rownames(res4_female), ensg.12k_female$V4),]
poss_dup<-ensg.12kb_female$V6=="IDS"
table(poss_dup)
#FALSE  TRUE 
#12813     2 
# IDS is duplicated in both males and females
ensg.dup<-ensg.12kb_female[poss_dup,]
# Look at expression patterns of the dups: 
summary(as.numeric(res4_female["ENSG00000010404.13",]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.0    53.5    78.0   114.9   118.0  1088.0  
summary(as.numeric(res4_female["ENSG00000241489.3",]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.2161  0.0000 15.0000 

# This is true for both females and males.
# Since the 2nd is lowly expressed, remove the second.. ENSG00000241489.3
ensg.to.rm_female<-ensg.12kb_female$V4=="ENSG00000241489.3"
ensg.to.rm_male<-ensg.12kb_male$V4=="ENSG00000241489.3"
ensg.12kc_female<-ensg.12kb_female[!ensg.to.rm_female,]
ensg.12kc_male<-ensg.12kb_male[!ensg.to.rm_male,]

res.to.rm_female<-rownames(res4_female)=="ENSG00000241489.3"
res.to.rm_male<-rownames(res4_male)=="ENSG00000241489.3"
res5_female<-res4_female[!res.to.rm_female,]
res5_male<-res4_male[!res.to.rm_male,]

# Let's double check again that the order is the same: 
ensg.12kd_female<-ensg.12kc_female[match(rownames(res5_female), ensg.12kc_female$V4),]
ensg.12kd_male<-ensg.12kc_male[match(rownames(res5_male), ensg.12kc_male$V4),]

res6_female<-res5_female
res6_male<-res5_male
rownames(res6_female)<-ensg.12kd_female$V6
rownames(res6_male)<-ensg.12kd_male$V6

setwd("/Users/kmagnaye/Documents/X_immunity_project/Processing_in_R")

write.table(res6_female, "GeneCounts_Female_GeneNames_noQC_12kgenes_110917.txt", sep="\t", quote=FALSE)
write.table(res6_male, "GeneCounts_Male_GeneNames_noQC_12kgenes_110917.txt", sep="\t", quote=FALSE)

# all expressed genes in the males are in the females and there are 22 female only genes

female_only <- (intersect(rownames(res6_female), rownames(res6_male))[!union(rownames(res6_female), rownames(res6_male)) %in% intersect(rownames(res6_female), rownames(res6_male))])
female_only
#  [1] "PRKX"    "FIGF"    "ZRSR2"   "DMD"     "LANCL3"  "SRPX"    "EFHC2"
# [8] "RBM10"   "SPIN3"   "OPHN1"   "CXorf65" "OGT"     "CITED1"  "BRWD3"  
# [15] "CENPI"   "RNF128"  "RGAG1"   "ZNF75D"  "CD40LG"  "MAMLD1"  "RENBP"
# [22] "EMD"  

######## 8. Split up files by treatment for downstream processing: ===============

lps_female <- pheno_7M[females & pheno_7M$treat=="LPS",]
null_female <- pheno_7M[females & pheno_7M$treat=="null",]
cd_female <- pheno_7M[females & pheno_7M$treat=="CD",]

lps_male <- pheno_7M[!females & pheno_7M$treat=="LPS",]
null_male <- pheno_7M[!females & pheno_7M$treat=="null",]
cd_male <- pheno_7M[!females & pheno_7M$treat=="CD",]

lps.gc_female <- res6_female[,colnames(res6_female) %in% lps_female[,2]]
null.gc_female <- res6_female[,colnames(res6_female) %in% null_female[,2]]
cd.gc_female <- res6_female[,colnames(res6_female) %in% cd_female[,2]]

lps.gc_male <- res6_male[,colnames(res6_male) %in% lps_male[,2]]
null.gc_male <- res6_male[,colnames(res6_male) %in% null_male[,2]]
cd.gc_male <- res6_male[,colnames(res6_male) %in% cd_male[,2]]

lps_female_match <- lps_female[lps_female[,2] %in% colnames(lps.gc_female),]
null_female_match <- null_female[null_female[,2] %in% colnames(null.gc_female),]
cd_female_match <- cd_female[cd_female[,2] %in% colnames(cd.gc_female),]

lps_male_match <- lps_male[lps_male[,2] %in% colnames(lps.gc_male),]
null_male_match <- null_male[null_male[,2] %in% colnames(null.gc_male),]
cd_male_match <- cd_male[cd_male[,2] %in% colnames(cd.gc_male),]

# Combine into two sets of treatments: 
nl.pheno_female<-rbind.data.frame(null_female_match, lps_female_match) # 70 + 77
nt.pheno_female<-rbind.data.frame(null_female_match, cd_female_match) # 70 + 67
nl.gc_female<-cbind.data.frame(null.gc_female, lps.gc_female) # 65 + 69
nt.gc_female<-cbind.data.frame(null.gc_female, cd.gc_female) # 65 + 64

nl.pheno_male<-rbind.data.frame(null_male_match, lps_male_match) # 54 + 54
nt.pheno_male<-rbind.data.frame(null_male_match, cd_male_match) # 54 + 54
nl.gc_male<-cbind.data.frame(null.gc_male, lps.gc_male) # 41 + 38
nt.gc_male<-cbind.data.frame(null.gc_male, cd.gc_male) # 41 + 40

setwd("/Users/kmagnaye/Documents/X_immunity_project/Processing_in_R/PreQC_GeneCount_pheno_files")

# female only files

write.table(null_female_match, "Female_Pheno_Rready_null_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(lps_female_match, "Female_Pheno_Rready_lps_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cd_female_match, "Female_Pheno_Rready_cd_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(nl.pheno_female, "Female_Pheno_Rready_nullLPS_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(nt.pheno_female, "Female_Pheno_Rready_nullCD_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(null.gc_female, "Female_RawGeneCounts_Null_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)
write.table(lps.gc_female, "Female_RawGeneCounts_LPS_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)
write.table(cd.gc_female, "Female_RawGeneCounts_CD_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)

write.table(nl.gc_female, "Female_RawGeneCounts_NullLPS_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)
write.table(nt.gc_female, "Female_RawGeneCounts_NullCD_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)

# Male only files

write.table(null_male_match, "Male_Pheno_Rready_null_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(lps_male_match, "Male_Pheno_Rready_lps_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(cd_male_match, "Male_Pheno_Rready_cd_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(nl.pheno_male, "Male_Pheno_Rready_nullLPS_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(nt.pheno_male, "Male_Pheno_Rready_nullCD_7M_110917.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(null.gc_male, "Male_RawGeneCounts_Null_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)
write.table(lps.gc_male, "Male_RawGeneCounts_LPS_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)
write.table(cd.gc_male, "Male_RawGeneCounts_CD_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)

write.table(nl.gc_male, "Male_RawGeneCounts_NullLPS_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)
write.table(nt.gc_male, "Male_RawGeneCounts_NullCD_RunsCombined_GeneNames_noQC_110917.txt", sep="\t", quote=FALSE)





