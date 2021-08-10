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
    """Gets coordinate of CA atoms in a chain and the indices which are None"""
    #list_of_conformation_coordinate_matrices
    #pymol indices are 1 --> 298, however 
    #this list will be 0 --> 297 for the sake of future programmability
    cmd.deselect()
    coords = list()
    for i in range(1,chain_length+1):
        coords.append(cmd.get_coords(f"/{chain}///{i}/CA", 1))
    return coords

def multiprocessing_get_coords(chains_list, chain_length=298):
    p = Pool()
    result = p.map(get_coords,chains_list)
    print("results retrieved")
    with open("coords_tup_res.var", "wb") as infile1:
        pickle.dump(result,infile1)
    with open("coords_tup_res.txt", "w") as infile2:
        infile2.write(str(result))
def pca_all_coordinates(all_chain_coords,chains_list):
    """Takes a list of chain coordinate matrices, add all the points into one massive matrix"""
    all_coords_list= list()
    for matrix in all_chain_coords:
        all_coords_list+= matrix
    all_coords_array = np.matrix(all_coords_list)
    ################ make sure to only sum every i of 298.
    avg = np.sum(all_coords_array, axis=0)/ 20 #change this into number of not none inputs of the coordinates


print(f"__name__: {__name__}")
main()