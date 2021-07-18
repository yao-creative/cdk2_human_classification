from pymol import cmd
import os
import pickle
from best_align import align_choice

def main():
    print(f"main running")
    with open("annotated.var","rb") as annotated_var:
        annotated_dict = pickle.load(annotated_var)
    with open("annotated_chains.var","rb") as annotated_chains_var:
        annotated_chains_dict= pickle.load(annotated_chains_var)
    choice = align_choice()
    print(f"choice: {choice}")
    #print(f"annotated_dict: {annotated_dict}")
    parent_code = choice[0][:-2]
    cmd.load(f"PDB_files/{parent_code.lower()}.pdb")
    cmd.select(f"/{parent_code}//{choice[0][-1]}")
    cmd.extract(f"{choice[0]}","sele")
    for chain in annotated_dict[parent_code].chains:
        if chain != choice[0][-1]:
            cmd.select(f"/{parent_code}//{chain}")
            cmd.extract(f"{parent_code}_{chain}","sele")
            cmd.align(f"{parent_code}_{chain}",choice[0])
            cmd.deselect()
            cmd.delete(parent_code)
    iterate_align(choice,parent_code,annotated_dict)
def iterate_align(choice,choice_parent,annotated_dict):
    for code in annotated_dict:
        conformation = annotated_dict[code]
        #print(f"conformation: {conformation}")
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
print("hi1")
main()       
print("hi")




