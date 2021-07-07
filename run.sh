#!/bin/bash

cd PDB_files
ls *.pdb
run=#?
cd ../
if [ $run -eq None ]; then
    echo "No pdb files, downloading pdb files"
    curl https://www.uniprot.org/uniprot/P63208.txt > human_cdk2_uniprot.txt
    python3 pdb_codes.py
    python3 pdb_get.py
    cd PDB_files
    ./batch_download.sh -f searchlist.txt -p
    python3 uncompress.py
    cd ../
fi
python3 annotate.py
open annotated.txt

