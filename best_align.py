import pickle
#although it would've been nice to add a class pdb dictionary and add this as a method, I realized too late
#and it would be too much of a hassle to rework all of the structure of the code
def align_choice(annotated_dict):
    res_list = list()
    for code in annotated_dict:
        for chain in annotated_dict[code].chains:
            print(f"align choice: {code}")
            conformation = annotated_dict[code]
            #print(f"code: {code} resolution: {conformation.res} atps: {len(conformation.atps)}")
            if code.lower() == "4eoj_a":
                print(f"4eoj_a: {repr(conformation)}")
            try:
                if float(conformation.res[:-1]) < 2.0 and chain in conformation.atps:
                    res_list.append((code, chain, conformation.res))
            except:
                continue
    print(f"res_list: {res_list}")
    choice = min(res_list, key = lambda tup: tup[2])
    return (f"{choice[0]}_{choice[1]}", choice[2])


#print(f"align_choice: {align_choice()}")
#print(f"align choice parent: {align_choice()[0][:-2]}")