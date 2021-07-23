from pymol import cmd
import os
import pickle
from best_align import align_choice
import numpy as np
from alignment_indices import collusion_list
def main():
    print(f"started")
    #load files
    with open("annotated.var","rb") as annotated_var:
        annotated_dict = pickle.load(annotated_var)
    print("loaded 1")
    annotated_var.close()
    #pick suitable choice to align the rest of the conformations to, choice will be a CDK2 chain of a conformation
    print(f"done loading chains!")
    choice = align_choice(annotated_dict)
    print(f"choice: {choice}")
    #load and align the chains of the chosen focus
    parent_code = choice[0][:-2]
    cmd.load(f"PDB_files/{parent_code.lower()}.pdb")
    cmd.color("sand", f"/{parent_code}") 
    cmd.select(f"/{parent_code}//{choice[0][-1]}")
    cmd.extract(f"{choice[0]}","sele")
    for chain in annotated_dict[parent_code].chains:
        if chain != choice[0][-1]:
            cmd.select(f"/{parent_code}//{chain}")
            cmd.extract(f"{parent_code}_{chain}","sele")
            cmd.align(f"{parent_code}_{chain}",choice[0])
            cmd.deselect()
            cmd.delete(parent_code)
    #iterate and align the rest of the conformations
    iterate_align(choice,parent_code,annotated_dict)
    remove_het()
    cmd.deselect()
    cmd.save("all_aligned.pse", "all_aligned")
    print(f"alignment state saved")
    cmd.deselect()
    #rmsd_matrix(annotated_dict)
    return

def iterate_align(choice,choice_parent,annotated_dict, res_threshold=(False,None)):
    """Takes the choice of alignment code, the parent of the choice, 
    which is already loaded and annotated dictionary of codes. option res_threshold
    to not load certain files with lower resolution (bool: True/False, float: min resolution)
    Returns list of skipped conformations"""
    skipped_list= list()
    for code in annotated_dict:
        conformation = annotated_dict[code]
        #print(f"conformation: {conformation}")
        if res_threshold[0]:
            if res_threshold[1] <= conformation.res:
                skipped_list.append(code)
                continue #conformation resolution too low, skip.
        print(f"loading code: {code}")
        if code == choice_parent: #exclude the chain part of the string
            print(f"found choice!")
            continue
        try:
            cmd.load(f"PDB_files/{code.lower()}.pdb") 
            cmd.color("sand", f"/{code}") 
            for chain in conformation.chains:
                cmd.select(f"/{code}//{chain}")
                cmd.extract(f"{code}_{chain}","sele")
                cmd.align(f"{code}_{chain}",choice[0])
                cmd.deselect()
            cmd.delete(f"{code}")
            cmd.remove(f"resn hoh")
        except:
            pass
    return skipped_list
 
def rmsd_matrix(chains_list):
    n = len(chains_list)
    order = list()
    matrix = np.zeros(shape=(n,n))
    cmd.alter("all", "segi=''")
    cmd.alter("all", "chain=''")
    for i,chain_id in enumerate(chains_list):
        for j,other_id in enumerate(chains_list):
            ############################
            #modify matrix entry
            matrix[i][j]= cmd.rms_cur(f"/{chain_id}////CA",f"/{other_id}////CA",matchmaker=4)
            print(f"chain: {chain_id} other: {other_id} rms success {matrix[i][j]}")  
    print(f"matrix: {matrix}")
    with open("matrix.var","wb") as infile1:
        pickle.dump(matrix,infile1)
        infile1.close()
    with open("matrix.txt","w") as infile2:
        for row in matrix:
            infile2.write(str(row))
        infile2.close()
        
            

        #assuming it goest through the dictionary in the same order
def remove_het():
    with open("chains_list.var","rb") as chains_list_var:
        chains_list = pickle.load(chains_list_var)
    for chain_id in chains_list:
        try:
            cmd.select(f"/{chain_id}///298:/")      
            cmd.remove("sele")
        except:
            pass
    cmd.remove("resn ATP")
    cmd.remove("resn MG")







print(f"started")
with open("chains_list.var","rb") as chains_list_var:
    chains_list= pickle.load(chains_list_var)
cmd.load("all_aligned.pse")
print(f"loaded all align!")
rmsd_matrix(chains_list) 
main()
print("hi")




