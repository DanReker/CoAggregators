#!/bin/bash
###
# Required input: a tab-delimited file such that first column conatins name of molecule, and second column contains SMILES string
# USAGE: 
#	1. chmod u+x run_obabel.sh
#	2. ./run_obabel.sh name_smiles.txt

#This code will run :
# 1. openbabel (create mol2 file with conformation & partial charges)
# 2. antechamber (create random residue name, change atom type into gaff format)
# 3. and parmchk2 (add additional atom/bond parameters)
###

while read name smiles; do
# give random residue name to molecules of interest in order to avoid ovrlap with existing residues
id=$(LC_CTYPE=C tr -dc A-Z < /dev/urandom | head -c 2 | xargs)
resname=Z$id

# First search & create optimal conformation systematically, if the process goes over 2hrs, just turn to weighted search
if timeout 7200 obabel -:"$smiles $name" -omol2 --conformer --gen3d --systematic -h --ff GAFF --partialcharge eem -O $name.mol2; then
  echo "$name - Systematic!"
else
  obabel -:"$smiles $name" -omol2 --conformer --gen3d --weighted -h --ff GAFF --partialcharge eem -O $name.mol2
  echo "$name - Weighted!"
fi

# Assure that atom types are recognized as gaff atom type; residue name is corrected; Atom names become unique
antechamber -i $name.mol2 -fi mol2 -o $name.mol2 -fo mol2 -rn $resname

# Add default value of parameters if certain bonds/atoms/angles are anomalies.
parmchk2 -i $name.mol2 -f mol2 -o $name.frcmod 

# if anything fails, create a buffer file that contains such information
if [ ! -s "$name.frcmod" ]; then
echo $name >> FAILED.txt
fi
done < $1

