import pickle
import math
import list_compare

def main(choice, choice_parent, method= "cl", dump= True):
    with open("annotated.var","rb") as annotated_var:
        annotated_dict = pickle.load(annotated_var)
    annotated_var.close()
    alignment_dict = get_alignment(choice, choice_parent,annotated_dict,method,dump)
    list_compare.main()

#################################################################

#aligment method 1 Needleman-Wunsch

def nw_alignment(seq1,seq2,cache=None):
    if cache is None: cache = (0,list())


#################################################################

#aligment method 2 collusion list

def collusion_list(a,b, cache= None):
    """given two lists A and B as arguments,
     returns a list of pairs of indices [(i1,j1),…,
     (ik,jk)] representing their longest matching 
     subsequence"""
    #the issue about this is that potentially there's a sequence
    #KVDAADEG and KVADE, however the A in the final AD of the second sequence matches 
    #first A however in reality is the second A? and then the indices are messed up and will be filtered out.
    if cache == None:
        cache = {}
    if (len(a), len(b)) in cache:
        return cache[len(a), len(b)]
    if len(a)<=0 or len(b)<=0:
        return []
    if a[-1]== b[-1]:
        #print("Equality!")
        cache[len(a), len(b)] = collusion_list(a[0:-1],b[0:-1],cache) + [(len(a)-1,len(b)-1)]
        return cache[len(a),len(b)]
    else: 
        #print("finding max")
        temp1= collusion_list(a[0:-1],b,cache)
        temp2= collusion_list(a,b[0:-1],cache)
        if len(temp1) > len(temp2):
            cache[len(a), len(b)] = temp1
        else: 
            cache[len(a), len(b)] = temp2
        return cache[len(a), len(b)]

#################################################################

#get the well aligned indices
def get_alignment(choice,choice_parent,annotated_dict,method="nw", dump = True):
    """takes the chosen sequence and dictionary of other conformations returns a list of alignment indices
    (choice,choice_parent,annotated_dict,method="nw", dump = True)"""
    print(f"getting alignment!")
    if method != "nw" and method != "cl":
        print(f"invalid method!")
        return
    choice_chain_residues = annotated_dict[choice_parent].residues[choice[-1]]
    print(f"choice_chain: {choice_chain_residues}")
    alignment_dict = dict()
    #first go through and find individual alignments
    for code in annotated_dict: #420
        conformation = annotated_dict[code]
        print(f"code: {code}")
        for chain in conformation.chains: #max size 2
            chain_residues= conformation.residues[chain]
            if f"{code}_{chain}" == choice:
                alignment_dict[choice] = [[i for i in range(len(chain_residues))],[i for i in range(len(chain_residues))]]
                continue #the rms_cur will be 0 of course
            else: #run dynamic algorithm to align the residue lists
                i = 0
                if method == "nw":
                    alignment_dict[f"{code}_{chain}"]= nw_alignment(choice_chain_residues,chain_residues)
                elif method == "cl":
                    tup_list = collusion_list(choice_chain_residues,chain_residues)
                    l1=list()
                    l2=list()
                    for tup in tup_list: #max size 298
                        l1.append(tup[0])
                        l2.append(tup[1])
                    alignment_dict[f"{code}_{chain}"]= [l1,l2]
    print(f"alignment_dict created!")
    #print(f"alignment_dict: {alignment_dict}")
    if dump:
        with open("alignment_dict.var","wb") as alignment_dict_var:
            pickle.dump(alignment_dict,alignment_dict_var)
        alignment_dict_var.close()
    print("alignment_dict dumped!")
    #print(f"dumped! alignment indices starting")
    #return alignment_dict

    
if __name__ == "__main__":
    
    main("4EOJ_A","4EOJ")

#print(collusion_list(["a","b", "c", "d", "d"], ["f","a","b","e","c"]))

#annotated_var.close()




