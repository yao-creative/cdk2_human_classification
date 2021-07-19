import os
import pickle
from Bio import PDB
from Bio import Align
#Non functionining program, ignore for the momet
def get_structure(annotated_dict):
    parser = PDB.PDBParser(PERMISSIVE= True, QUIET= True)

    for item in annotated_dict:
        item =item.lower()
        structure = parser.get_structure(item, f"{item}.pdb")
        
        for model in structure:
            residues = dict()
            for chain in model:
                if chain.id in annotated_dict[item.upper()].chains:
                    chain_res= list()#keys: chains, values: residue in list
                    for residue in chain:
                        if residue.resname != "HOH":
                            chain_res.append(residue.resname)
                    residues[chain.id] = chain_res
                #print('>some_header\n',''.join(seq))
        """models =list(structure.get_models())
        print(f"structure.get_models(): {models}")
        print(f"models[0]: {models[0]}")
        residues = [list(chain.get_residues) for chain in [list(structure.get_models().get_chains())]]"""
        atoms = [list(residue) for residue in residues]
        conformation =annotated_dict[item.upper()]
        conformation.structure= structure
        conformation.residues= residues
        conformation.atoms = atoms
        annotated_dict[item.upper()]= conformation
        print(f"{item} success!\nResidue: {residues}")
    
