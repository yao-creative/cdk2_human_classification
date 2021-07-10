from pymol import cmd
import os
import pickle
import sys

print(f"str(sys.argv): {str(sys.argv)}")

groups =dict()
closed_inactive= open("/Users/yao/Desktop/dkp/Work/internship1_bioinfo/structures/closed_active.var", "rb")
opened_active = open("/Users/yao/Desktop/dkp/Work/internship1_bioinfo/structures/opened_active.var", "rb")
opened_inactive = open("/Users/yao/Desktop/dkp/Work/internship1_bioinfo/structures/opened_inactive.var", "rb")

groups["opened_active"]=pickle.load(opened_active)
groups["closed_inactive"]=pickle.load(closed_inactive)
groups["opened_inactive"]=pickle.load(opened_inactive)

os.chdir("../PDB_files")
print(f"current: {os.getcwd()}")

#cmd.load(f"{groups['opened_active'][0]}.pdb")
group_list = [code.lower() for code in groups[str(sys.argv[2])]]
for i,code in enumerate(group_list):
    print(f"code: {code}")
    try:
        cmd.load(f"{code}.pdb")
    except:
        pass
    if i !=0: #align all the other structures to the first one
        try:
            cmd.align(code, group_list[0])
        except:
            pass

        print(f"aligned: {code} to {group_list[0]}")

print(f"group_list: {group_list}")
#os.system("pause")
#input("press any key to continue")
