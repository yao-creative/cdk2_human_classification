import json


#CLASSES
#---------------------------------------------------
#
class Pdb:
    def __init__(self, code, type, res, chains, length, group=None, structure=None, atoms =None):
        self.code = code
        self.type = type
        self.res = res
        self.chains = chains
        self.length =length
        self.file =None
        self.group = group
        #using store bio python structure feature as one of it's attributes
        self.structure = structure
        self.atoms = atoms
    def __repr__(self, structure=False):
        if structure:
            return f"code: {self.code}; type: {self.type} res: {self.res} chains: {self.chains} length: {self.length} group: {self.group}\nStructure: {self.structure}\natoms: {self.atoms}"
        else:
            return f"code: {self.code}; type: {self.type} res: {self.res} chains: {self.chains} length: {self.length} group: {self.group}"






#FUNCTIONS
#---------------------------------------------------
#
def nxt_itm(line,idx):
    """function which takes a line, and index the start of a series of spaces
    and returns the string right after the series of spaces before the next series of spaces"""
    k=idx
    while line[k] == " " and k<len(line):
        k+=1
    j = k
    while line[j] != " " and k<len(line):
        j+=1
    return line[k:j],j

def dumper(obj):
    try:
        return obj.toJSON()
    except:
        return obj.__dict__