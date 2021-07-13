from pymol import cmd
import os
import pickle
import sys
sys.path.insert(1, "../")
import helper
print(f"str(sys.argv): {str(sys.argv)}")


groups =dict()
closed_inactive= open("/Users/yao/Desktop/dkp/Work/internship1_bioinfo/structures/closed_active.var", "rb")
opened_active = open("/Users/yao/Desktop/dkp/Work/internship1_bioinfo/structures/opened_active.var", "rb")
opened_inactive = open("/Users/yao/Desktop/dkp/Work/internship1_bioinfo/structures/opened_inactive.var", "rb")

groups["opened_active"]=pickle.load(opened_active)
groups["closed_inactive"]=pickle.load(closed_inactive)
groups["opened_inactive"]=pickle.load(opened_inactive)
os.chdir("..")
with open("annotated.var","rb") as annotated_var:
    annotated_dict = pickle.load(annotated_var)
#print(f"annotated_dict: {annotated_dict}")
#print(f"annotated_dict[1e9h]: {annotated_dict['1E9H'].chains}")
os.chdir("PDB_files")
print(f"current: {os.getcwd()}")

#cmd.load(f"{groups['opened_active'][0]}.pdb")
group_list = [code.lower() for code in groups[str(sys.argv[2])]]
for i,code in enumerate(group_list):
    print(f"code: {code}")
    
    try:
        cmd.load(f"{code}.pdb") 
        cmd.color("sand", f"/{code}") 
        chain =annotated_dict[code.upper()].chains[0]
        cmd.remove(f"{code} and (not Chain {chain})")
        cmd.remove(f"resn hoh")
        #cmd.color("deep salmon",f"{code} & Chain A & index 33-44")
        cmd.select(f"/{code}//{chain}/33:44/CA")
        cmd.color("marine",f"sele")
        cmd.deselect()
        cmd.select(f"/{code}//{chain}/9:17/CA")
        cmd.color("tv_red","sele")
    except:
        pass
    if i==0:
        focus = code
        print(f"focus: {focus}")
    if i !=0: #align all the other structures to the first one
        try:
            cmd.align(code, focus)
        except:
            pass

        print(f"aligned: {code} to {focus}")
print(f"group_list: {group_list}")
#os.system("pause")
#input("press any key to continue")
