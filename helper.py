import pickle
import os

#CLASSES
#---------------------------------------------------
#
class Pdb:
    def __init__(self, code, type, res, chains, length, group=None, structure=None, atoms =None, residues=None, rms_cur=None, atps=[]):
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
        self.residues = residues
        self.rms_cur = rms_cur
        self.atps = atps
    def __repr__(self, show=None,):
        if show == "rms":
            return f"rms_cur: {self.rms_cur}"
        if show == "structure":
            return f"code: {self.code}; type: {self.type} res: {self.res} chains: {self.chains} length: {self.length} group: {self.group} atps: {self.atps}\nStructure: {self.structure}\nresidues: {self.residues}"
        else:
            return f"code: {self.code}; type: {self.type} res: {self.res} chains: {self.chains} length: {self.length} group: {self.group} atps: {self.atps}"








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