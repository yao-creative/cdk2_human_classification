from pymol import cmd
import os
import pickle
import numpy as np
#from alignment_indices import collusion_list
from multiprocessing import Pool
import time
import sys
#from pml_script_all import align_all,get_chains_list#,comparision_mat
from helper import bin_search
from sklearn.decomposition import PCA
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

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

    #cmd.load("all_aligned.pse")
    print(f"all aligned loaded")
    if "coords_res.var" not in lsdir:
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
    p.join()
    p.close()
    print("results retrieved")
    with open("coords_res.var", "wb") as infile1:
        pickle.dump(result,infile1)
    with open("coords_res.txt", "w") as infile2:
        infile2.write(str(result))
def pca_all_coordinates(all_chain_coords,chains_list):
    """Takes a list of chain coordinate matrices, add all the points into one massive matrix"""
    all_coords_list= list()
    for matrix in all_chain_coords:
        all_coords_list+= matrix
    all_coords_array = np.matrix(all_coords_list)
    ################ make sure to only sum every i of 298.
    avg = np.sum(all_coords_array, axis=0)/ 20 #change this into number of not none inputs of the coordinates

##### Normally there'd be some variation of the gradient ascent algorithm and threshold to find high variation segments, 
# for the purpose of time I'll impute values gotten from experimentation.

def pick_seg(coords):
    pass
def reduce_seg(coords,seg_idx_list,chains_list,threshold=10):
    """Takes in the coordinates and picks out important segments
    seg_idx_list list of tuples of start and end 0 indexed
    list of chains, and threshold if a segment has more than 5 nones then the whole chain is thrown away"""
    seg_list = list() #list of conformation list of segment matrices which contains actual residue coordinates 
    #list(list(matrix()))
    for matrix in coords:
        temp = list() #list of segment matrix new matrix
        for seg in seg_idx_list:
            temp.append(matrix[seg[0]:seg[1]+1])
        seg_list.append(temp)
    removed_indices = list()
    for i,conformation in enumerate(seg_list):
        stop = False
        for seg_mat in conformation:
            #print(f"seg_mat: {seg_mat}")
            if sum(1 for i in seg_mat if i is None) >= threshold:
                removed_indices.append(i)
                break
    
    print(f"removed_indices: {removed_indices}")
    new_seg_list = list()
    reduced_chains_list= list()
    for i,item in enumerate(seg_list):
        #print(f"i: {i} ")
        if bin_search(removed_indices,0, len(removed_indices)-1,i):
            continue
        else:
            new_seg_list.append(item)
            reduced_chains_list.append(chains_list[i])
    #print(f"removed indices: {removed_indices} new_seg_list: {new_seg_list} reduced_chains_list: {reduced_chains_list}")
    return new_seg_list,reduced_chains_list
def impute(matrix):
    imp = IterativeImputer()
    return imp.fit_transform(matrix)

def coordinate_impute(seg_list):
    """Impute missing values
    Parameters: seg_list: as list of conformations which are represented by a list of segment coordinate matrices
    list(list(matrix()))"""
    conformation_by_seg = list() #list of segments within each segment are the coordinates 
    #for all of conformations corresponding to that segment
    for i in range(len(seg_list[0])):
        conformation_by_seg.append([conformation[i] for conformation in seg_list])
    print(f"length conformation by seg[0]: {len(conformation_by_seg[0])}")
    conformation_by_seq_by_index = list()
    show = True
    for segment in conformation_by_seg:
        temp = [list() for _ in range(len(segment[0]))]
        for conformation in segment:
            for i, residue in enumerate(conformation):
                if residue is None:
                    residue = np.array([np.nan, np.nan, np.nan])
                else:
                    residue = residue[0]
                
                temp[i].append(residue)
        for i in range(len(temp)):
            #print(f"temp[i]: {temp[i]}")
            #print(f"shape temp[i]: {np.shape(temp[i])}")
            temp[i] = np.matrix(temp[i])
        p = Pool()
        result =  p.map(impute,temp)
        p.close()
        p.join()
        conformation_by_seq_by_index.append(result)
    with open("imputed_seg_coord.var","wb") as infile:
        pickle.dump(conformation_by_seq_by_index,infile)
    with open("imputed_seg_coord.txt","w") as infile2:
        infile2.write(str(conformation_by_seq_by_index))
    return conformation_by_seq_by_index
    #print(f"conformation by seq by index: {conformation_by_seq_by_index}")


def seg_val(conformation_by_seq_by_index):
    """takes seg_list and reduces it into one feature,
    coords list of matrices of coordinates
    Parameters: list(list(np.matrix()))
    A list of segment lists which contain conformation coordinates matrix for that segment"""

    #Unpack: 
    
    #print(conformation_by_seg)
    
    # all_coords_mat = list()
    # for matrix in seg_list:
    #     all_coords_mat+= matrix
    ##Map all of the coordinates single values.
    pca = PCA(n_components=1)
    #all_coords_reduced = pca.fit_transform(all_coords_mat)

    
with open("coords_res.var", "rb") as infile3:
    coords = pickle.load(infile3)
with open("chains_list.var","rb") as chains_list_var:
    chains_list= pickle.load(chains_list_var)
#seg_list,reduced_chains_list = reduce_seg(coords,[(32,43),(149,158)],chains_list)
# with open("reduced_chains_list.var", "wb") as reduced_chains_list_var:
#     pickle.dump(reduced_chains_list,reduced_chains_list_var)
# with open("seg_list.var","wb") as seg_list_var:
#     pickle.dump(seg_list,seg_list_var)
#     seg_list_var.close()
def run():
    with open("seg_list.var","rb") as seg_list_var:
        seg_list= pickle.load(seg_list_var)
        seg_list_var.close()
    with open("reduced_chains_list.var", "rb") as reduced_chains_list_var:
        reduced_chains_list= pickle.load(reduced_chains_list_var)
        reduced_chains_list_var.close()

    #coordinate_impute(seg_list)
    with open("imputed_seg_coord.var","rb") as infile:
        conformation_by_seq_by_index= pickle.load(infile)
        infile.close()
    
run()
#seg_val(seg_list)
#print(f"seg_list: {seg_list}")


print(f"__name__: {__name__}")
#main()