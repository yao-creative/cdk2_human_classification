import pickle
import os
outfile = open("pdbs2.txt", "w")
os.chdir("..")
with open("pdbs.var", "rb") as f:
    pdb_dict = pickle.load(f)


for pdb in pdb_dict:
    outfile.write(pdb + " ")

outfile.close()
f.close()
