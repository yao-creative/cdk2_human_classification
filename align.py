import os
import pickle
from Bio import PDB
from Bio import Align

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'} 
parser = PDB.PDBParser(PERMISSIVE= True, QUIET= True)
with open("pdbs.var", "rb") as f:
    pdb_dict = pickle.load(f)
os.chdir("PDB_files")

for item in pdb_dict:
    item =item.lower()
    structure = parser.get_structure(item, f"{item}.pdb")
    
    for model in structure:
        for chain in model:
            residues = []
            for residue in chain:
                residues.append(d3to1[residue.resname])
            print('>some_header\n',''.join(seq))
    """models =list(structure.get_models())
    print(f"structure.get_models(): {models}")
    print(f"models[0]: {models[0]}")
    residues = [list(chain.get_residues) for chain in [list(structure.get_models().get_chains())]]"""
    atoms = [list(residue) for residue in residues]
    conformation =pdb_dict[item]
    conformation.structure= structure
    conformation.residues= residues
    conformation.atoms = atoms
    pdb_dict[item]= conformation
    
