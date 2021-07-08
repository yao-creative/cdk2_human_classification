import os
import pickle
from Bio import PDB
from Bio import Align

parser = PDB.PDBParser(PERMISSIVE= True, QUIET= True)
with open("pdbs.var", "rb") as f:
    pdb_dict = pickle.load(f)
os.chdir("PDB_files")

for item in pdb_dict:
    item =item.lower()
    pdb_dict[item].structure = parser.get_structure(item, f"{item}.pdb")
    pdb_dict[item].atoms = [list(residue) for residue in [list(chain.get_residues) for chain in [list(model.get_chains) for model in list(pdb_dict[item].structure.get_models())]]]
    
