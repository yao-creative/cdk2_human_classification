import numpy as np
import pickle
chains_list2 = list()
with open("distmat.dat", "r") as infile:
    all_list = list()
    indices = list()
    for line in infile:
        l = line.split()
        all_list.append(l)
        indices.append(int(l[0]))
    n = max(indices)
    print(f"n: {n}")
    matrix = np.zeros(shape=(n+1,n+1))
    for l in all_list:
        matrix[int(l[0]),int(l[1])] = float(l[-1])
        if f"{l[2][:-1]}_{l[2][-1]}" not in chains_list2:
            chains_list2.append(f"{l[2][:-1]}_{l[2][-1]}")
    
    print(f"chains_list2: {chains_list2}")
    infile.close()
with open("chains_list2.var", "wb") as infile2:
    pickle.dump(chains_list2,infile2)
    infile2.close()
matrix = matrix + matrix.T
print(f"matrix: {matrix}")
with open("dihedral_mat.var", "wb") as infile3:
    pickle.dump(matrix,infile3)
    infile3.close()
