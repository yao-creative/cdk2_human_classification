# CDK2 Classification
This project is based on the classif_cdk2_conf.pdf file in reference papers where I wrote scripts to first annotate the groups of 
Open inactive, open active, closed inactive using the pdb files of uniprot P24941 CDK2, based on the presence of TPO and cyclin structures.
Then using RMSD and high variance regions method I cluster the proteins and check for the statistics. The full detail is explained in the paper: Exploring_Intuitive_Approaches_to_Protein_Conformation_Clustering_Using_Regions_of_High_Structural_Variance.pdf. 
## Note:
Anywhere it says coordinate-based, refers to the T-SNE Based map on the paper

## Description:
Goes to the internet retrieves uniprot codes of Human CDK2 conformations, runs a script to obtain the pdb files and model them on pymol while aligning all of them to one of the conformation then outputs RMSD matrix of the conformation vs other conformations.
## Requirements:
Python3, Pymol, Sklearn, scipy, pickle, matplotlib
## Instructions:
1) Download the uniprot file from uniprot which has the list of pdbfiles you want to get.
2) go to terminal and type ```bash run.sh```
if you want to reset and test another protein, run reset.sh
1) the annotations script has to be modified 
2) The name of the groups with the pymol lists have to be modified
3) The pml_script_all has to be modified to find specific RMSD values or requirements based on experiments. (only the cmd.rmsd() part) and also the files to save to.



## Python File descriptions:

### helper.py: 
contains pdb class and other useful functions

### annotate.py:
Annotates the structures pdb code by pdb code and stores it's properties
### annotate_chains.py:
Annotates the structures chain by chain code and stores it's properties. Run after annotate.py.

### reset.sh: 
Clears all the txt,var,pdb files.

### run.sh:
Runs all of the files in the correct order (make sure for different protein annotations have to be modified)

### pdb_codes.py:
 Program which iterates through the text version of the uniprot page and creates a list of pdb files

### pdb_get.py:
 takes the list of codes and extracts the pdb codes as objects and outputs out: a string of codes separated by commas and spaces and pdb_list a dictionary of pdb code objects dumping the dictionary in pdbs.var and prints a non annotated file pdbs.txt.

### look_through.py:
 Check through if every pdb file has specific type property

### uncompress.py:
 uncompresses all gz files in pdb_files

### pml_script_all.py:
pymol script which gets the best fitting structure based on criteria, iterates and aligns the rest of the structures to it, then saves the file all_aligned.pse. Afterwards it uses multiprocessing to calculate the RMSD distances based on criterias. Has to be modified if only finding RMSD for certain subsequences. Make sure to check the writing out commands and the cmd.rmsd(structure1, structure2) to see what it's doing.
### all_aligned.pse:
A pymol file where all of the proteins are aligned on top of the choice of structures. 
### pml_script_coord.py:
Has run at least after pml_script_all.py or at least after pml_script_all.py iterates and aligns all of the structures and saves a file all_aligned.pse. Retrieves the post-alignment coordinates of the structures after pymol is done.
### parse_dat.py:
Reads the PSS distance data and puts it into a usable matrix for clustering

### coordinate_cluster.ipynb:
Python notebook using the t-SNE based map clustering method, experimental results

### dihedral_clustering.ipynb:
Python notebook using the dihedral-based clustering method after being parsed from the PSS run, experimental results.

### rmsd_cluster.ipynb:
Python notebook using the rmsd matrix clustering method gotten from the pml_script_all.py, experimental results.

### find_hv_regions.ipynb:
Python notebook using the residue to next residue vector drawing method to find high variance substructures, experimental results.

### pymol_list.py:
Creates usable lists each group using chain-wise annotations from annotation_chains in the structures folder 
### Notes:
1) pml_script_coord is specific to CDK2 since it specifically takes into account coordinates of the CA from 1 to 297
2) all of the matrix*.var or .txt files are just experimental results fround my pml_script_coord or pml_script.