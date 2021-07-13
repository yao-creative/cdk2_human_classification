# CDK2 Classification

## Instructions:
### (if pdb files already downloaded skip first 5 steps)

### 1) Go to uniprot to download the text version of the human cdk2 proteins entry list.
### 2) Run pdb_list.py which gets CDK2 structures from uniprot txt
### 3) Run pdb_get.py which creates a list of obj PDB and also PDB codes in a txt usable format for the batch_download.sh
### 4) go into PDB_files and run ```./batch_download.sh -f searchlist.txt -p``` on the terminal to get the gzip files.
### 5) To unzip the files run uncompress.py

### 6) To annotate the files, run annotate.py and it'll out put annotated.txt an annotated version of the files


## Python File descriptions:

### helper.py: contains pdb class and other useful functions

### pdb_codes.py: Program which iterates through the text version of the uniprot page and creates a list of pdb files

### pdb_get.py: takes the list of codes and extracts the pdb codes as objects and outputs out: a string of codes separated by commas and spaces and pdb_list a dictionary of pdb code objects dumping the dictionary in pdbs.var and prints a non annotated file pdbs.txt.

### look_through.py: Check through if every pdb file has specific type property

### uncompress.py: uncompresses all gz files in pdb_files

## Requirements:
### Python3 
### Pymol