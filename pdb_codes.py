"""Program which iterates through the text version of the uniprot page and creates a list of pdb files"""

k = 2
#list of pdb objects
outfile = open("list.txt", "w")
#list of pdb codes separated by commas

with open("human_cdk2_uniprot.txt", 'r') as infile:
    read_pdb = True
    start = False
    for line in infile:
        #print(f"line: {line}")
        if line[:2] == "DR" and read_pdb:
            if line[5:9] == "PDB;":
                start = True
                outfile.write(line)
            elif start:
                break
outfile.close()
infile.close()