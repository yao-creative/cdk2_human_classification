from pymol import cmd
import os
import pickle
import numpy as np
#from alignment_indices import collusion_list
from multiprocessing import Pool
import time
import sys
from pml_script_all import align_all,get_chains_list#,comparision_mat

def main():
    #Note if we want to redo alignment, makes sure that the file all_aligned.pse is removed
    print(f"started")
    lsdir = os.listdir(".")
    #print(f"lsdir: {lsdir}")
    if "all_aligned.pse" not in lsdir:
        align_all()
    if "chains_list.var" not in lsdir:
        get_chains_list()    
    with open("chains_list.var","rb") as chains_list_var:
        chains_list= pickle.load(chains_list_var)

    cmd.load("all_aligned.pse")
    print(f"all aligned loaded")
    multiprocessing_get_coords(chains_list)
    print(f"program finished")

def get_coords(chain,chain_length=298):
    """Gets coordinate of CA atoms in a chain"""
    #list_of_conformation_coordinate_matrices
    coord_mat = list()
    existing_coords = list() #pymol indices are 1 --> 298, however 
    #this list will be 0 --> 297 for the sake of future programmability
    for i in range(1,chain_length+1): #This part is specific 
        try:
            cmd.deselect()
            cmd.select(f"/{chain}///{i}/CA")
            coord_mat.append(cmd.get_coords('sele', 1))
            existing_coords.append(i-1)
        except:
            continue
    return (existing_coords,coord_mat)

def multiprocessing_get_coords(chains_list, chain_length=298):
    p = Pool()
    result = p.map(get_coords,chains_list[0:2])

    print("results retrieved")
    with open("coords_tup_res.var", "wb") as infile1:
        pickle.dump(result,infile1)
    with open("coords_tup_res.txt", "w") as infile2:
        infile2.write(str(result))
    
def pca_all_coordinates():
    pass

print(f"__name__: {__name__}")
main()