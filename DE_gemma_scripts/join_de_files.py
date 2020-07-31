# -*- coding: utf-8 -*-
"""
testing locations

folder_loc = "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/testing_data/"
gender = "females"
chr_type = "X"

"""

#import sys module
import sys


def join_gemma(folder_loc, gender, chr_type):
    
    #get file names
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(folder_loc) if isfile(join(folder_loc, f)) and gender in f and ("." + chr_type + ".") in f and 'fe' + gender not in f and "all" not in f]
    
    #Define output list
    output_list =[]
    #loop through files and merge lists
    for i in onlyfiles:
        
        #Define empty list to take values to append to output list and gene holder
        temp_list = []
        temp_gene = i[:]
        #Add gene name to start of temp_list
        temp_gene = temp_gene.replace(gender + ".", "")
        temp_gene = temp_gene.replace(".gemma.assoc.txt", "")
        temp_gene = temp_gene.replace("." + chr_type, "")
        temp_list.append(temp_gene)
        
        #Read in i file
        hold_list =[]
        with open(folder_loc + "/" + i) as results_i:
            for line in results_i.readlines():
                hold_list.append(line.strip().split("\t"))         
        del line
        
        #Add remaining info to temp_list
        temp_list.extend(hold_list[1])
        
        #Add temp_list to output_list
        output_list.append(temp_list)
    
    #Write to file
    write_path = folder_loc + "/" + gender + ".all." + chr_type + ".gemma.assoc.txt"
    import csv
    with open(write_path, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
    
        for i in output_list:
            writer.writerow(i)
        writeFile.close()
        
#call function with system arguments input from command shell
join_gemma(sys.argv[1], sys.argv[2], sys.argv[3])