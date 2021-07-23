import pickle
with open("pdbs.var", "rb") as infile:
    dictionary = pickle.load(infile)
out = list()
infile2= open("chains_list.txt", "w")
for code in dictionary:
    for chain in dictionary[code].chains:
        out.append(f"{code}_{chain}")
        infile2.write(f"{code}_{chain}\n") 
with open("chains_list.var", "wb") as infile3:
    pickle.dump(out,infile3)
infile.close()
infile2.close()
infile3.close()