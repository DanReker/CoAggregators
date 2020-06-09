###
# Usage: run_OpenMM.sh directory_name INFILE (same input file from step 2) 
#   	- Creates a directory that contains the simulation files and copy of the simulation parameter python codes.
#	
# This code will create :
#	1. simulation trajectory file in pdb format
#	2. create log file that contains simulation statistics,
#	3. (optional) time_log file that records the length of computational time the simulation took to create
#
# To change simulation parameters, modify line 36-38.
###
  
#!/bin/bash
dir=$1
if [ ! -e $dir ]; then
mkdir $dir
fi
cp ./run_OpenMM.sh $dir/

# get filename of mixed box
file=$2
simfile_name=''
# get molecule identities from pair file
n_col=$(awk '{print NF}' ${file} | sort -nu | tail -n 1)
n_row=$(cat ${file} | wc -l)
for i in $(seq 1 ${n_row}); do
  for j in $(seq 1 $(( n_col / 2 ))); do
    mol=$(awk -v row=$i -v col=$j 'FNR==(row) {print $(col)}' ${file})
    simfile_name=${simfile_name}_${mol}
  done
  echo ${simfile_name:1}
  simfile_name=''
done > simfiles


# Modify this part to run OpenMM with wanted parameters
while IFS='' read -r name; do
  python OpenMM.py -top ${name}.prmtop -crd ${name}.inpcrd -o ${dir}/openmm_${name}.pdb -i
done < simfiles
