###
# Usage: ./run_OpenMM.sh simulation_directory ../data/pair_composition.tsv 
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


# Modify this part to run OpenMM with wanted parameters. More details can be shown by "python OpenMM.py -h"
# --savechk             If activated, save a checkpoint state at every 10 intervals
# --loadchk             If activated, load a checkpoint state with .chk extension of the output filename
# -c {0,1,2,3}          Bond constraints; 0 = None (default); 1 = HBonds ; 2 = AllBonds ; 3 =HAngles
# -p, --periodic        If activated, runs the simulation in periodic box with PME method used for long range interaction (default = NoCutoff)
# -i, --implicit        If activated, runs the simulation in implicit water solvent (default = vacuum, unless explicitly solvated)
# --solvent_type {1,2,3,4,5} 1 = HCT (default); 2 = OBC1 ; 3 = OBC2 ; 4 = GBn ; 5 = GBn2
# -T TEMPERATURE        Set simulation temperature (default = 300K)
# --timestep TIMESTEP   Set simulation time step in units of picosecond (default = 0.002 picosecond)
# --interval INTERVAL   Set interval of saved frames in the unit of picosecond (default = 10 ps)
# -l SIMULATION_LENGTH  Set duration of simulation time in units of nanosecond (default = 20ns)
# --timeit              If activated, creates time_log file
while IFS='' read -r name; do
  python OpenMM.py -top ${name}.prmtop -crd ${name}.inpcrd -o ${dir}/openmm_${name}.pdb -i
done < simfiles