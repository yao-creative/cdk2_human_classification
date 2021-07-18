"""after annotate.py, split each of the multichain conformation into their chains"""
import pickle
import helper
import copy
annotated_chains_txt  = open("annotated_chains.txt", "w")

with open("annotated.var","rb") as annotated_var:
    annotated_dict = pickle.load(annotated_var)
annotated_chains_dict = dict()
for code in annotated_dict:
    conformation = annotated_dict[code]
    for chain in conformation.chains:
        #print(f"conformation chains: {conformation.chains}")
        conformation_chain = copy.deepcopy(conformation)
        conformation_chain.chains = [chain]

        if chain in conformation_chain.atps:
            conformation_chain.atps = [chain]
        else:
            conformation_chain.atps = []
        conformation_chain.code = f"{code}_{chain}"
        annotated_chains_txt.write(repr(conformation_chain) + "\n")
        annotated_chains_dict[f"{code}_{chain}"] = conformation_chain

with open("annotated_chains.var", "wb") as f:
    pickle.dump(annotated_chains_dict, f)

annotated_chains_txt.close()
f.close()
annotated_var.close()    



