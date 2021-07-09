import pickle
import os
import sys
sys.path.insert(1, "../")
import helper
pymol_list = open("pymol_list.txt", "w")
opened_active = open("opened_active.txt", "w")
closed_inactive = open("closed_active.txt", "w")
opened_inactive = open("opened_inactive.txt", "w")

os.chdir("..")
with open("annotated.var", "rb") as f:
    pdb_dict = pickle.load(f)



for pdb in pdb_dict:
    pdb =pdb_dict[pdb]
    pymol_list.write(pdb.code + " ")
    if "open" in pdb.group:
        if "active" in pdb.group:
            opened_active.write(pdb.code + " ")
        else:
            opened_inactive.write(pdb.code + " ")
    else:
        closed_inactive.write(pdb.code + " ")
pymol_list.close()
f.close()
