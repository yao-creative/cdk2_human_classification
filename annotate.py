import os
from helper import nxt_itm as nxt_itm
import pickle 
from get_structure import get_structure

annotated_txt = open("annotated.txt", "w")

#retrieve dictionary of pdbs using pickle 
with open("pdbs.var", "rb") as f:
    pdb_dict = pickle.load(f)
os.chdir("PDB_files")
for item in os.listdir("."):
    #creates list of chains with TPO occurence 
    tpo_list = list()
    atp_list = list()
    if item.endswith("pdb"):
        #print(f"item name: {item}")
        with open(item, "r") as infile:
            #first classify by in active using tyrosine 160 phosphorylated or not.
            start = False
            opened = False #based on cyclin presence
            for line in infile:
                line= line.strip()
                if line[:6] == "COMPND":
                    temp =line.lower()
                    if temp.find("cyclin")>=0 and temp.find("cyclin-dependent")<0:
                        #make sure it's cyclin and not a cdk
                        opened = True
                if line[:4] == "HET ": #check heterogenous structure for TPO 160
                    start = True
                    molecule,i = nxt_itm(line,4)
                    if molecule == "TPO":
                        #if it's TPO check the chain identity
                        chain_id,j = nxt_itm(line,i)
                        seq_num,k=nxt_itm(line,j)
                        #makesure to only record tpos at position 160 of the chain
                        if seq_num == "160":
                            tpo_list.append(chain_id)
                    if molecule == "ATP":
                        chain_id,j = nxt_itm(line,i)
                        atp_list.append(chain_id[0])
                elif start:
                    #end after reading the formulas
                    break
            
            pdb =  pdb_dict.get(item[:-4].upper())
            pdb.atps = atp_list
            pdb.tpo_list = tpo_list
            #go through atps and the one 
            #no tpos
            if len(tpo_list) == 0:
                pdb.group = ["inactive"]
                if opened:
                    pdb.group.append("open")
                else: pdb.group.append("closed")
            #check if tpos are part of the right chain
            for chain_id in tpo_list:
                if chain_id in pdb.chains:
                    pdb.group = ["active", "open"]
                    break
                else: 
                    pdb.group = ["inactive"]
                    if opened:
                        pdb.group.append("open")
                    else: pdb.group.append("closed")
                
            pdb_dict[item[:-4].upper()] = pdb
            #print(f"pdb: {pdb}")

get_structure(pdb_dict)
for pdb in pdb_dict:
    annotated_txt.write(repr(pdb_dict[pdb]) + "\n")
os.chdir("..")
with open("annotated.var", "wb") as f:
    pickle.dump(pdb_dict,f)  #some how after getting the structures the dump has become heavy duty
f.close() 
annotated_txt.close()