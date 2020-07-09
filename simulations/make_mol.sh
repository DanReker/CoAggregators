#!/bin/bash
###
# Required input: a tab-delimited file such that first column contains name of molecule, and second column contains SMILES string
# USAGE: 
#	1. chmod u+x run_obabel.sh
#	2. ./make_mol.sh ../data/selected_drug_smiles.txt
#
# This code will execute:
# 1. openbabel (create mol2 file with conformation & partial charges)
# 2. Ambertools - antechamber (assign unique residue name, change atom type into gaff format)
# 3. Ambertools - parmchk2 (check and add additional atom/bond parameters)
###

while read name smiles; do
# Assign random residue name to molecules of interest in order to avoid ovrlapped parameters with existing residues.
id=$(LC_CTYPE=C tr -dc A-Z < /dev/urandom | head -c 2 | xargs) # this is only a hack, and there could be a more systematic way of generating unique IDs for each molecule. 
resname=Z$id

# Systematically search & create for the optimal conformation, but if the process goes over 2hrs (7200s), abort and run weighted search.
# Timeout duration can be adapted from the original 7200 sec
if timeout 7200 obabel -:"$smiles $name" -omol2 --conformer --gen3d --systematic -h --ff GAFF --partialcharge eem -O ${name}.mol2; then 
  echo "$name - Systematic!"
else
  obabel -:"$smiles $name" -omol2 --conformer --gen3d --weighted -h --ff GAFF --partialcharge eem -O ${name}.mol2
  echo "$name - Weighted!"
fi

# Assure that atom types are recognized as gaff atom type; residue name is corrected.
antechamber -i $name.mol2 -fi mol2 -o $name.mol2 -fo mol2 -rn $resname

# Add default value of parameters if certain bonds/atoms/angles are anomalies.
parmchk2 -i $name.mol2 -f mol2 -o $name.frcmod 

# if a molecules fails to generate, append in a buffer file.
if [ ! -s "$name.frcmod" ]; then
echo $name >> FAILED.txt
fi
done < $1