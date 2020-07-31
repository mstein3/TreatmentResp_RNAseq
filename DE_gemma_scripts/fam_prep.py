"""
Testing locations
expr_path = "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/testing_data/fam_files/females.ARSD.auto.txt"
fam_path = "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/testing_data/fam_files/Hutterite.females.rigged.auto.fam"
"""

#import sys module
import sys

def fam_prep(expr_path, fam_path, write_path):
    
    #Parse the pheno file into a list
    expr_raw = []
    with open(expr_path) as expr_file:
        for line in expr_file.readlines():
            expr_raw.append(line.strip().split("\t"))
    del line

    #Parse the fam file into a list
    fam_raw = []
    with open(fam_path) as fam_file:
        for line in fam_file.readlines():
            fam_raw.append(line.strip().split(" "))
    del line
    
    for i in range(len(fam_raw)):
        fam_raw[i][5] = expr_raw[i][2]

    #Write fam to file
    import csv
    with open(write_path, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = ' ')
        for i in fam_raw:
            writer.writerow(i)
        writeFile.close()


        
#call function with system arguments input from command shell
fam_prep(sys.argv[1], sys.argv[2], sys.argv[3])