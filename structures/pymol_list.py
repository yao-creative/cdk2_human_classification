import pickle
import os
import sys
sys.path.insert(1, "../")
import helper
pymol_list = open("pymol_list.txt", "w")
opened_active = open("opened_active.txt", "w")
closed_inactive = open("closed_active.txt", "w")
opened_inactive = open("opened_inactive.txt", "w")
opened_active_var = open("opened_active.var", "wb")
closed_inactive_var = open("closed_active.var", "wb")
opened_inactive_var = open("opened_inactive.var", "wb")

os.chdir("..")
with open("annotated.var", "rb") as f:
    pdb_dict = pickle.load(f)

opened_active_list = list()
closed_inactive_list = list()
opened_inactive_list = list()
for pdb in pdb_dict:
    pdb =pdb_dict[pdb]
    pymol_list.write(pdb.code + " ")
    if "open" in pdb.group:
        if "active" in pdb.group:
            opened_active.write(pdb.code + " ")
            opened_active_list.append(pdb.code)
        else:
            opened_inactive.write(pdb.code + " ")
            opened_inactive_list.append(pdb.code)
    else:
        closed_inactive.write(pdb.code + " ")
        closed_inactive_list.append(pdb.code)

pickle.dump(opened_active_list, opened_active_var)
pickle.dump(opened_inactive_list, opened_inactive_var)
pickle.dump(closed_inactive_list, closed_inactive_var)


pymol_list.close()
opened_active.close()
closed_inactive.close()
opened_inactive.close()
opened_active_var.close()
closed_inactive_var.close()
opened_inactive_var.close()
pymol_list.close()
f.close()
