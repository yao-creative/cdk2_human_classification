import pickle


with open("annotated.var", "rb") as infile1:
    annotated_dict = pickle.load(infile1)
infile1.close()
for code in annotated_dict:
    annotated_dict[code].structure= None
with open("annotated.var", "wb") as outfile1:
    pickle.dump(annotated_dict,outfile1)
outfile1.close()


with open("annotated_chains.var", "rb") as infile2:
    annotated_chains_dict = pickle.load(infile2)
infile2.close()
for code in annotated_chains_dict:
    annotated_chains_dict[code].structure= None
with open("annotated_chains.var", "wb") as outfile2:
    pickle.dump(annotated_chains_dict,outfile2)
outfile2.close()
