import os
os.chdir("PDB_files")
focus = ""
count = 0
for item in os.listdir("."):
    #creates list of chains with TPO occurence 
    if item.endswith(".pdb"):
        #print(f"item name: {item}")
        stop = False
        with open(item, "r") as infile:
            for line in infile: