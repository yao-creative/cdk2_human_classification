#!/bin/bash

cd PDB_files
ls *.pdb
get_pdbs=$?
#validity is 0 if there are pdb
echo "validity $get_pdbs"
if [ $get_pdbs != 0 ] ; then
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