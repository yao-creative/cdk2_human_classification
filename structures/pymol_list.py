import pickle
import os
import sys
sys.path.insert(1, "../")
import helper
def main():
    print("started pymol_list.py")
    pymol_list = open("pymol_list.txt", "w")
    opened_active = open("opened_active.txt", "w")
    closed_inactive = open("closed_inactive.txt", "w")
    opened_inactive = open("opened_inactive.txt", "w")
    opened_active_var = open("opened_active.var", "wb")
    closed_inactive_var = open("closed_inactive.var", "wb")
    opened_inactive_var = open("opened_inactive.var", "wb")

    os.chdir("..")
    with open("annotated_chains.var", "rb") as f:
        chains_dict = pickle.load(f)

    opened_active_list = list()
    closed_inactive_list = list()
    opened_inactive_list = list()
    for chains_code in chains_dict:

        conformation_chain =chains_dict[chains_code]
        pymol_list.write(f"{conformation_chain.code}" + " ")
        if "open" in conformation_chain.group:
            if "active" in conformation_chain.group:
                opened_active.write(conformation_chain.code + " ")
                opened_active_list.append(conformation_chain.code)
            else:
                opened_inactive.write(conformation_chain.code + " ")
                opened_inactive_list.append(conformation_chain.code)
        else:
            closed_inactive.write(conformation_chain.code + " ")
            closed_inactive_list.append(conformation_chain.code)

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
    
if __name__ == "__main__":
    main()