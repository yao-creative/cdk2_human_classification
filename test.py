import pickle
with open("matrix_AB.var", "rb") as infile1:
    matrix = pickle.load(infile1)
print(f"matrix[0] len: {len(matrix[0])}")