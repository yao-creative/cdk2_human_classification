# CDK2 Classification


## Description:
Goes to the internet retrieves uniprot codes of Human CDK2 conformations, runs a script to obtain the pdb files and model them on pymol while aligning all of them to one of the conformation then outputs RMS matrix of the conformation vs other conformations.
## Requirements:
Python3, Pymol
## Instructions:
### go to terminal and type ```bash run.sh```



## Python File descriptions:

### helper.py: contains pdb class and other useful functions

### pdb_codes.py: Program which iterates through the text version of the uniprot page and creates a list of pdb files

### pdb_get.py: takes the list of codes and extracts the pdb codes as objects and outputs out: a string of codes separated by commas and spaces and pdb_list a dictionary of pdb code objects dumping the dictionary in pdbs.var and prints a non annotated file pdbs.txt.

### look_through.py: Check through if every pdb file has specific type property

### uncompress.py: uncompresses all gz files in pdb_files

## Requirements:
### Python3 
### Pymol