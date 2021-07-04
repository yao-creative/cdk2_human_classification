import os
from helper import Pdb  as Pdb
from helper import dumper as dumper
import json
import pickle

#create empty pdb list of pdb objects
pdbs_txt = open("pdbs.txt", "w")

os.chdir("PDB_files")
ltxt = open("searchlist.txt", "w")

os.chdir("..")
#sort into pdb objects
def unravel_chain(chain_str):
    chains = list()
    for char in chain_str:
        if char == "/":
            pass
        else: chains.append(char)
    return chains
def make_pdb_obj_list():
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
            res = line[idx[2]+1: idx[2]]
            eq = line.index("=")
            chain_str = line[idx[3]+1: eq]
            chains = unravel_chain(chain_str)
            length = line[eq+3:len(line)-2]
            curr =Pdb(code,type,res,chains,length)
            pdb_dict[code] =  curr
            out += code.lower()+", "
            pdbs_txt.write(repr(curr) + "\n")
    
    infile.close()
    return out, pdb_dict

out, pdb_dict =make_pdb_obj_list()
#create list seperated by commas of pdb codes
ltxt.write(out)



#Failed attempt to use pickle to dump my protein object variables.
print(f"pdb_dict: {pdb_dict}")
with open("pdbs.var", "wb") as f:
    pickle.dump(pdb_dict,f)
    #print(f"dumped: {pdb_dict}")
#f.close()


pdbs_txt.close()
ltxt.close()

