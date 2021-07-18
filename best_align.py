import pickle
#although it would've been nice to add a class pdb dictionary and add this as a method, I realized too late
#and it would be too much of a hassle to rework all of the structure of the code
def align_choice():
    with open("annotated_chains.var","rb") as annotated_chains_var:
         annotated_chains_dict= pickle.load(annotated_chains_var)

    res_list = list()
    for code in annotated_chains_dict:
        conformation = annotated_chains_dict[code]
        #print(f"code: {code} resolution: {conformation.res} atps: {len(conformation.atps)}")
        if code.lower() == "4eoj_a":
            print(f"4eoj_a: {repr(conformation)}")
        try:
            if float(conformation.res[:-1]) < 2.0 and len(conformation.atps) >0:
                res_list.append((code, conformation.res))
        except:
            continue
    print(f"res_list: {res_list}")
    choice = min(res_list, key = lambda tup: tup[1])
    return choice


#print(f"align_choice: {align_choice()}")
#print(f"align choice parent: {align_choice()[0][:-2]}")