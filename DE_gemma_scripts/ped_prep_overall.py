"""
Testing locations

male_pheno_loc = "C:/Users/Janel/Documents/sex-specific/processed_data/males.Pheno.Rready.nullLPS.7M.031520.txt"
female_pheno_loc = "C:/Users/Janel/Documents/sex-specific/processed_data/females.Pheno.Rready.nullLPS.7M.031520.txt"
grm_id_path = "C:/Users/Janel/Documents/sex-specific/filters/overall.iids.txt"

"""

#import sys module
import sys

def ped_prep_grant(female_pheno_loc, male_pheno_loc, grm_id_path, ped_path, map_path):
    
    #Parse the pheno files into a list
    pheno_raw = []
    with open(female_pheno_loc) as pheno_file:
        for line in pheno_file.readlines():
            pheno_raw.append(line.strip().split("\t"))
    del line
    del pheno_raw[:1]
    temp = []
    with open(male_pheno_loc) as pheno_file:
        for line in pheno_file.readlines():
            temp.append(line.strip().split("\t"))
    del line
    del temp[:1]
    pheno_raw.extend(temp)
    del temp
    
    #Open and list of grant iids
    id_raw = []
    with open(grm_id_path) as id_file:
        for line in id_file.readlines():
            id_raw.append(line.strip().split("\t"))
    #Cycle through id_raw and remove first entry in each list
    for i in range(len(id_raw)):
        id_raw[i] = int(id_raw[i][1])
    #Sort id_raw, though order should already be correct
    id_raw.sort()
    del i
    
    
    #Need to cycle through pheno_raw and convert extract ords to ints
    for entry in pheno_raw:
        entry[1] = int(entry[1])
    del entry
    #Sort pheno_raw
    pheno_raw.sort(key = lambda x:x[1])


    #cycle through entries in ref_list and ped_raw to create ped_final
    #Also orders correctly to match pheno file
    ped_final = []
    for i in range(len(pheno_raw)):
        if pheno_raw[i][1] in id_raw:
            #create temp_list to hold values    
            temp_list = []
            #Add first five column values to temp list
            temp_list.append('Hutterites')
            temp_list.append(pheno_raw[i][1])
            temp_list.extend(pheno_raw[i][76:80])
            #Add rigged columns for treatment
            if pheno_raw[i][4] == 'null':
                temp_list.extend(['A', 'A'])
            elif pheno_raw[i][4] == 'LPS':
                temp_list.extend(['T', 'T'])
            ped_final.append(temp_list)
    #clean-up
    del i
    del pheno_raw
    del temp_list
    del id_raw
    
    #sort ped_final to match GRM order
    ped_final.sort(key = lambda x:x[1])

    map_final = [[23, 9999, 0, 9999]]    
        
    #Write ped to file
    import csv
    with open(ped_path, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = ' ')
        for i in ped_final:
            writer.writerow(i)
        writeFile.close()

    #Write ped to file
    import csv
    with open(map_path, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
        for i in map_final:
            writer.writerow(i)
        writeFile.close()
        
#call function with system arguments input from command shell
ped_prep_grant(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])