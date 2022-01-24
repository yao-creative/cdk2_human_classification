import os
#Check through if every pdb file has specific type property
os.chdir("PDB_files")
focus = "MOLECULE: CYCLIN-"
count = 0
for item in os.listdir("."):
    #creates list of chains with TPO occurence 
    
    if item.endswith(".pdb"):
        #print(f"item name: {item}")
        found = False
        with open(item, "r") as infile:
            for line in infile:
                #print(f"line: {line}")
                res = line.find(focus)
                if res>=0:
                    found = True
                    #print(f"item: {item} res: {res}")
                    break
        
        if found: 
            #if the focus wasn't found in the file 
            #notify not consistent and end program.
            print(f"{item} found")
            count+=1
        else:
            print("{item} not found")
        
print(f"count: {count}")
            
        