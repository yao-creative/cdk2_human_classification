#!/bin/bash

cd PDB_files
ls *.pdb
get_pdbs=$?
#validity is 0 if there are pdb
echo "validity $get_pdbs"
if [ $get_pdbs != 0 ] ; then
    cd ../
    echo "No pdb files, downloading pdb files"
    curl https://www.uniprot.org/uniprot/P24941.txt > human_cdk2_uniprot.txt
    python3 pdb_codes.py
    python3 pdb_get.py
    cd PDB_files
    ./batch_download.sh -f searchlist.txt -p
    python3 uncompress.py
    
fi
cd ../
echo "annotating"
python3 annotate.py
open annotated.txt


#render the groups on pymol
cd structures
ls *.txt
split_groups=$?
if [ $split_groups != 0 ] ; then
    python3 pymol_list.py
fi

pymol pml_script.py opened_active
pymol pml_script.py closed_inactive
pymol pml_script.py opened_inactive