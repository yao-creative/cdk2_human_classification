from pymol import cmd
import os
import pickle
from best_align import align_choice

def main():
    print(f"main running")
    #load files
    with open("annotated.var","rb") as annotated_var:
        annotated_dict = pickle.load(annotated_var)
    
    #pick suitable choice to align the rest of the conformations to, choice will be a CDK2 chain of a conformation
    choice = align_choice()
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




    
def rmsd_matrix():
    pass
print("hi1")
main()       
print("hi")




