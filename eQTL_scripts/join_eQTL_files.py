# -*- coding: utf-8 -*-
"""
testing locations

folder_loc = "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/testing_data/eQTL_testing"
chr_type = "auto"

"""

#import sys and math modules
import sys
import math
import random

'''This function is an interior function that takes an in progress list of gene results, reads 
    in a list of results for a new condition, matches rows and adds them to the correct place'''
def read_and_match(gene_list, folder_loc, gene, chr_type, gender, treatment):
    #Read in file and sort it by rs column
    raw_list =[]
    with open(folder_loc + "/" + "eQTL." + gender + "." + gene + "." + chr_type + "." + treatment + ".gemma.assoc.txt") as results_i:
        for line in results_i.readlines():
            raw_list.append(line.strip().split("\t"))         
    del line
    raw_list = raw_list[1:]
    raw_list.sort(key = lambda x:x[1])
    
    #Create return_list to actually return data
    #This allows for discarding of non-matches
    return_list = []
    #Create starting position tracking variable for while loop
    start_pos = 0
    #Loop through gene_list 
    for i in range(len(gene_list)):
        #Reset match to 0
        match = 0
        #Reset current_pos to start_pos
        current_pos = start_pos
        #Check match condition and only then exit inner loop and jump to next map_variant
        while match == 0:
            #Check match
            if gene_list[i][2] == raw_list[current_pos][1]:
                #Set match to 1 and update start_pos
                match = 1
                start_pos = current_pos
                #Fill temp_list
                temp_list = gene_list[i][:]
                temp_list.extend(raw_list[current_pos][7:9])
                temp_list.append(float(raw_list[current_pos][14]))
                #Send temp_list to x_variant_processed
                return_list.append(temp_list)
            #Check no match condition where map_variant has already been passed
            #by x_variant_id loop
            elif gene_list[i][2] < raw_list[current_pos][1]:
                #Set match to 1 to escape loop, but no match is actually found
                match = 1
            #Need to add a check to make sure not at the end of the raw_list
            #If this is the case it means that the new file has fewer snps than the old file,
            #and all of the snps in the new file have been used
            elif current_pos == len(raw_list) - 1:
                #Introducing a match condition will cause all future snps in gene_list to be excluded
                match = 1
            #No match and not yet passed location
            else:
                current_pos = current_pos + 1
    #Clean-up
    del match
    del current_pos
    del i
    del start_pos
    #Return Return_list
    return return_list
        

'''This function will do the following:
    1. Make a list of all genes tested in all four conditions
    2. Prep output list to receive info on total test count and beta's, SEs and p-values
        a. Also sort gene list to speed up future matching
    3. Cycle through genes and fill columns'''
def join_gene_lists(folder_loc, chr_type):
    
    #get file names
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(folder_loc) if isfile(join(folder_loc, f)) and ("." + chr_type + ".") in f and "combined" not in f]
        
    #Get complete list of genes
    genes = []
    for file in onlyfiles:
        #Get break positions of genes
        genes.append(file[file.find("males.") + 6 : file.find("males.") + 6 + file[file.find("males.") + 6:].find(".")])
    genes = list(dict.fromkeys(genes))
    
    #Create vectors of genders and treatment conditions to be useful later in loops
    genders = ["males", "females"]
    treatments = ["treatment", "no.treatment"]
    
    #Verify that gene was tested for eQTLs in all conditions
    verified_genes = []
    for gene in genes:
        #Create match variable that switches to off if a match isn't found
        match = True
        for gender in genders:
            for treatment in treatments:
                if not "eQTL." + gender + "." + gene + "." + chr_type + "." + treatment + ".gemma.assoc.txt" in onlyfiles:
                    match = False
        if match == True:
            verified_genes.append(gene)
    #Clean up
    del gene
    del gender
    del match
    del treatment
    del genes
    del file
    del onlyfiles
    
    
    '''Prep for loop'''
    #Define output list
    output_list =[]
    #Define total number of tests and make first line in output_list
    output_list.append(list())
    output_list[0] = [0]
    #Define headers of ouput_list
    output_list.append(['gene', 'chr', 'rs', 'ps', 'n_miss', 'allele1', 'allele0', 'af', 'ML_beta', 'ML_se', 'ML_p_score'
                         , 'MV_beta', 'MV_se', 'MV_p_score', 'FL_beta', 'FL_se', 'FL_p_score', 'FV_beta', 'FV_se', 'FV_p_score', 'Min_p_score'])
    #Sort verified gene list
    verified_genes.sort(key = lambda x:x[0])
    
    '''Create outer loop to cycle through genes in verified_genes'''
    for gene in verified_genes:
        #Define gene_list, which will hold results of individuals genes before sampling
        gene_list = []
        #Read in males_LPS file to fill first 10 columns of table
        raw_list =[]
        with open(folder_loc + "/" + "eQTL." + genders[0] + "." + gene + "." + chr_type + "." + treatments[0] + ".gemma.assoc.txt") as results_i:
            for line in results_i.readlines():
                raw_list.append(line.strip().split("\t"))         
        del line
        raw_list = raw_list[1:]
        #Add all males_LPs info to gene_list
        for entry in raw_list:
            if float(entry[6]) >= 0.05:
                temp_list = [gene]
                temp_list.extend(entry[0:9])
                temp_list.append(float(entry[14]))
                gene_list.append(temp_list)
                del temp_list
        #Clear junk
        del raw_list
        del entry
        #Verify that at least some of the variants had maf > 0.05
        if len(gene_list) > 0:
            #Sort gene list in preparation for function calls
            gene_list.sort(key = lambda x:x[2])
            #Function calls
            gene_list = read_and_match(gene_list, folder_loc, gene, chr_type, genders[0], treatments[1])
            gene_list = read_and_match(gene_list, folder_loc, gene, chr_type, genders[1], treatments[0])
            gene_list = read_and_match(gene_list, folder_loc, gene, chr_type, genders[1], treatments[1])
            #Cycle through gene list rows and add a minimum p-value to each row
            for entry in gene_list:
                entry.append(min(entry[10], entry[13], entry[16], entry[19]))
            del entry
            #Sort gene_list
            gene_list.sort(key = lambda x:x[20])
        
            #Add number of tests to sum in first row of output_list
            output_list[0][0] = output_list[0][0] + len(gene_list)
            #Append lead snp to output_list
            output_list.append(gene_list[0])
            #Get number of random snps to pull (i.e. 2% rounded down but with a minimum of 1)
            num_snps = int(max(math.floor(len(gene_list) * 0.02), 1))
            #Check if size of gene_list is one
            if len(gene_list) > 1:
                #Append a random sample of remaining snps to output_list
                output_list.extend(random.sample(gene_list[1:], num_snps))
            #Give out message to gague progress
            print(gene + " is complete.")
    
    #Clean-up
    del gene
    del gene_list
    del treatments
    del verified_genes
    del genders
    del num_snps
    
    #Write to file
    write_path = folder_loc + "/" + "combined." + chr_type + ".gemma.assoc.tsv"
    import csv
    with open(write_path, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
        for i in output_list:
            writer.writerow(i)
        writeFile.close()
    
#call function with system arguments input from command shell
join_gene_lists(sys.argv[1], sys.argv[2])
    