import os
from helper import Pdb  as Pdb
import json
import pickle

def main():
    pdbs_txt = open("pdbs.txt", "w")

    os.chdir("PDB_files")
    ltxt = open("searchlist.txt", "w")

    os.chdir("..")
    out, pdb_dict =make_pdb_obj_list(pdbs_txt)
    #create list seperated by commas of pdb codes
    ltxt.write(out)

    #Failed attempt to use pickle to dump my protein object variables.
    print(f"pdb_dict: {pdb_dict}")
    with open("pdbs.var", "wb") as f:
        pickle.dump(pdb_dict,f)
        #print(f"dumped: {pdb_dict}")
    #f.close()
    print(f"number of samples: {len(pdb_dict)}")


    pdbs_txt.close()
    ltxt.close()
#create empty pdb list of pdb objects

#sort into pdb objects
def unravel_chain(chain_str):
    chains = list()
    for char in chain_str:
        if char == "/":
            pass
        else: chains.append(char)
    return chains
def make_pdb_obj_list(pdbs_txt):
    """Function which takes the list of codes and extracts the pdb codes as objects and outputs 
    out: a string of codes separated by commas and spaces and pdb_list a list of pdb code objects"""
    pdb_dict = dict()
    out = ""
    
    with open("list.txt", "r") as infile:
        
        for line in infile:
            if line[0] == " ":
                continue
            line = line.replace(" ", "")
            idx = [i for i, x in enumerate(line) if x == ";"]
            #print(f"line: {line}\nidx: {idx}")
            code = line[idx[0]+1: idx[1]]
            type = line[idx[1]+1: idx[2]]
            res = line[idx[2]+1: idx[3]]
            eq = line.index("=")
            chain_str = line[idx[3]+1: eq]
            chains = unravel_chain(chain_str)
            length = line[eq+3:len(line)-2]
            len_str = ""
            for char in length:
                try:
                    len_str+=char
                    int(len_str) #make sure length is an integer
                except:
                    break
            curr =Pdb(code,type,res,chains,len_str)
            #pdb_dictionary with code as the key
            pdb_dict[code] =  curr
            out += code.lower()+", "
            pdbs_txt.write(repr(curr) + "\n")
    
    infile.close()
    return out, pdb_dict
if __name__ == "__main__":
    main()


