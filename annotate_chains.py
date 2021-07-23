"""after annotate.py, split each of the multichain conformation into their chains"""
import pickle
import helper
import json
import copy

def main():
    annotated_chains_txt  = open("annotated_chains.txt", "w")
    with open("annotated.var","rb") as annotated_var:
        annotated_dict = pickle.load(annotated_var)
    print(f"loaded!")
    annotated_chains_dict = dict()
    for code in annotated_dict:
        print(f"processing: {code}")
        conformation = annotated_dict[code]
        for chain in conformation.chains:
            #print(f"conformation chains: {conformation.chains}")
            conformation_chain = copy.deepcopy(conformation)
            conformation_chain.chains = [chain]

            if chain in conformation_chain.atps:
                conformation_chain.atps = [chain]
            else:
                conformation_chain.atps = []
            conformation_chain.residues = conformation.residues[chain] #might not be needed 
            #since already stored as dictionary in the annotated_dict
            conformation_chain.code = f"{code}_{chain}"
            annotated_chains_txt.write(repr(conformation_chain) + "\n")
            annotated_chains_dict[f"{code}_{chain}"] = conformation_chain

    #The following dump takes awhile as all the other codes too and might not be useful if we dump the residues
    with open("annotated_chains.var", "wb") as f:
        pickle.dump(annotated_chains_dict, f)

    annotated_chains_txt.close()
    f.close()
    annotated_var.close()    

if __name__ == "__main__":
    main()


