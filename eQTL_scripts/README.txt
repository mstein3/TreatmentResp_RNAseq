eQTL data temporarily stored in /scratch/mstein3/eQTL_Hutterites

Rough order of scripts: 

A. Before running eQTL: 
# Set up output directories 

B. Running eQTL 

run_short_GeneLists.sh 
# Submits a separate job for each gene list

RunGemma_ShGeneL_MN.pbs
RunGemma_ShGeneL_ML.pbs
RunGemma_ShGeneL_FN.pbs
RunGemma_ShGeneL_FL.pbs
# Submits the eQTL job for gemma for each dataset of interest (MN: Male untreated; ML: Male LPS-treated; FN: Female untreated; FL: Female LPS-treated)

Gemma_eQTL_perGene.R
# Runs gemma for a single gene, pulling from a matrix containing all genes. Genotypes in PLINK format. 

C. Getting top results by z-score 
Run_Assemble_top_results.pbs
# Submits job to calculate top z-scores

parse_gemma_topZscore.R
# Calculates z scores for each file, then calculates with SNP has the highest absolute value z-score for each gene (outputs files as *topZ.txt). 
# Finally, it looks at all of those *topZ.txt files, and combines them into a single file that has the top Z score for each gene within each data set. 

D. Calculate the top z-score SNP across datasets 
# Done in R. Generate a SNPlist that contains the top SNP (by z-score, across all datasets) for each Gene. 

E. Pull out the results from that curated SNP list from D. : 





run_short_GeneLists.sh 
Run_MaxBetaPull.pbs
MaxBetaPull.R

Assemble_MaxBetas.R 

