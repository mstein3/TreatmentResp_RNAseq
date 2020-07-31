"""
Testing locations

gene = "ACAP3"
gender = 'males'
gene_key = "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/filters/gene_locs_autosome.tsv"
map_file = "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/testing_data/Hutterite.males.process.auto.map"
write_loc = "C:/Users/Mitch V/Documents/UChicago/Research/X-linked_traits/testing_data/"

"""

#import sys module
import sys

def make_snv_list(gene, gender, gene_key, map_file, write_loc):
    #Parse the gene key into a list
    key_raw = []
    with open(gene_key) as key_file:
        for line in key_file.readlines():
            key_raw.append(line.strip().split("\t"))
    del line
    #Parse the map file into a list
    map_raw = []
    with open(map_file) as snv_file:
        for line in snv_file.readlines():
            map_raw.append(line.strip().split("\t"))
    del line

    #Create list to hold output
    output_list = []

    #Find starting position of gene of interest in gene key
    for entry in key_raw:
        if entry[0] == gene:
            #Recode non-numeric chromosomes
            if entry[1] == 'X':
                entry[1] = 23
            elif entry[1] == 'Y':
                entry[1] = 24
            elif entry[1] == 'M':
                entry[1] = 25
            gene_pos = int(entry[1]) * 10000000000 + int(entry[2])
            #cycle through map_raw to fill output_list if the gene is found in the key
            for entry2 in map_raw:
                #Recode non-numeric chromosomes
                if entry2[0] == 'X':
                    entry2[0] = 23
                elif entry2[0] == 'Y':
                    entry2[0] = 24
                elif entry2[0] == 'M':
                    entry2[0] = 25
                #Check 1Mb range around TSS
                if gene_pos - 1000000 <= int(entry2[0]) * 10000000000 + int(entry2[3]) and gene_pos + 1000000 >= int(entry2[0]) * 10000000000 + int(entry2[3]):
                    output_list.append([entry2[1]])
            del entry2
    del entry
    del key_raw
    

    #make output file path
    write_path = write_loc + gene + "." + gender +".snvs.txt"
    #Write pheno to file
    import csv
    with open(write_path, "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = ':')
        for i in output_list:
            writer.writerow(i)
        writeFile.close()

#call function with system arguments input from command shell
make_snv_list(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])