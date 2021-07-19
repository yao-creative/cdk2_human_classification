import pickle
def main():
    with open("alignment_dict.var","rb") as alignment_dict_var:
        alignment_dict = pickle.load(alignment_dict_var)
    out=get_alignment_indices(alignment_dict)
    with open("common_indices.txt","w") as infile1:
        infile1.write(str(out))
        infile1.close()
    with open("common_indices.var","wb") as infile2:
        pickle.dump(out,infile2)
        infile2.close()
##########################################################################

def bin_search(arr, l, h, x):
    # Check base case
    if h >= l:
        m = (h+ l) // 2
        if arr[m] == x:
            return True
        elif arr[m] > x:
            return bin_search(arr, l, m- 1, x)
        else:
            return bin_search(arr, m + 1, h, x)
    else:
        return False

def binary_list_elimination(values):
    """Takes in a list of lists and outputs a sublist where all the lists have in common
    n^2logn comparison speed"""
    #print(f"values length: {len(values)}")
    if len(values)== 1:
        return [values[0]]
    if len(values) == 2:
        out = []
        if type(values[1]) == int:
                print(f"values[1]: values[1]")
        for i in values[0]:
            if bin_search(values[1],0,len(values[1])-1,i):
                out.append(i)
        return [out]
    else:
        m = len(values)//2
        l1 = binary_list_elimination(values[:m])
        l2 = binary_list_elimination(values[m:])
        out = l1+ l2
        return binary_list_elimination(out)
def get_alignment_indices(alignment_dict):
    list_of_lists = list()
    for code in alignment_dict:
        if alignment_dict.get(code)[0] == 0:
            print(f"code: {code} 0 found!") 
        list_of_lists.append(alignment_dict.get(code)[0])
    print(f"list of lists: {list_of_lists} length: {len(list_of_lists)}")
    return binary_list_elimination(list_of_lists)
if __name__ == "__main__":
    main()
    