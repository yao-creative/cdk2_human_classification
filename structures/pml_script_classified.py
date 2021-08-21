from pymol import cmd
import os
import pickle
import sys
sys.path.insert(1, "../")
import helper


def main():
    print(f"main started")
    groups =dict()
    closed_inactive= open("closed_inactive.var", "rb")
    opened_active = open("opened_active.var", "rb")
    opened_inactive = open("opened_inactive.var", "rb")

    groups["opened_active"]=pickle.load(opened_active)
    groups["closed_inactive"]=pickle.load(closed_inactive)
    groups["opened_inactive"]=pickle.load(opened_inactive)
    os.chdir("..")
    with open("annotated.var","rb") as annotated_var:
        annotated_dict = pickle.load(annotated_var)
    #print(f"annotated_dict: {annotated_dict}")
    #print(f"annotated_dict[1e9h]: {annotated_dict['1E9H'].chains}")
    print(f"files loaded")
    os.chdir("PDB_files")
    #cmd.load(f"{groups['opened_active'][0]}.pdb")
    group_list = [code.lower() for code in groups[str(sys.argv[2])]]
    loaded = list()
    print(os.listdir("."))
    for i,code in enumerate(group_list):
        print(f"code: {code}")
        code = code[:-2]
        
        #print(f"code processed: {code}\n loaded: {loaded}")
        #print(f"loaded: {loaded}")
        if code in loaded:
            #print("continued")
            continue
        loaded.append(code)
        
        #try:
        cmd.load(f"{code}.pdb") 
        #print(f"loaded: {code}")
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
        cmd.deselect()
        cmd.select(f"/{code}//{chain}/150:159/CA")
        cmd.color("palegreen","sele")
        # except:
        #     pass
        if i==0: #we're expecting that each group_list is non-empty
            focus = code
            print(f"focus: {focus}")
        if i !=0: #align all the other structures to the first one
            try:
                cmd.align(code, focus)
            except:
                pass

            print(f"aligned: {code} to {focus}")
    cmd.save(f"{str(sys.argv[2])}.pse", "all_aligned")
    print(f"group_list: {group_list}")
    

main()
#os.system("pause")
#input("press any key to continue")
