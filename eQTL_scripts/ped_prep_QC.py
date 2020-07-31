"""
Testing locations

ped_loc = "C:/Users/Mitch V/Documents/UChicago/Research/"
pheno_loc = "C:/Users/Mitch V/Documents/UChicago/Research/Sex-specific/processed_data/males.Pheno.Rready.nullLPS.7M.031520.txt"
write_path = "C:/Users/Mitch V/Documents/UChicago/Research/"

"""

#import sys module
import sys

def ped_prep(ped_loc, pheno_loc, write_path):
    #Parse the ped file into a list
    ped_raw = []
    with open(ped_loc) as UTR_file:
        for line in UTR_file.readlines():
            ped_raw.append(line.strip().split(" "))
    del line
    #Parse the pheno file into a list
    pheno_raw = []
    with open(pheno_loc) as pheno_file:
        for line in pheno_file.readlines():
            pheno_raw.append(line.strip().split("\t"))
    del line

    #make ref list of ids to convert ped ids to pheno ids
    #Create list of extract ids from pheno file
    ref_list = []
    for entry in pheno_raw[1:]:
        ref_list.append([int(entry[2]), int(entry[1])])
    del entry

    #sort ped_final to match GRM order
    ref_list.sort(key = lambda x:x[1])
    
    #cycle through entries in ref_list and ped_raw to create ped_final
    #Also orders correctly to match pheno file
    ped_final = []
    temp_list = []
    for i in range(len(ref_list)):
        for j in range(len(ped_raw)):
            #Check if findivs match
            if ref_list[i][0] == int(ped_raw[j][1]):
                #create temp list to modify iid column and then append to ped_final
                temp_list = ped_raw[j][:]
                temp_list[1] = ref_list[i][1]
                ped_final.append(temp_list)
    #clean-up
    del i
    del j
    del pheno_raw
    del ref_list
    del temp_list
    del ped_raw
            
    #Write pheno to file
    import csv
    with open(write_path, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = ' ')
        for i in ped_final:
            writer.writerow(i)
        writeFile.close()

#call function with system arguments input from command shell
ped_prep(sys.argv[1], sys.argv[2], sys.argv[3])