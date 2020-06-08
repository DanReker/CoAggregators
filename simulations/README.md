# Co-Aggregators Simulations

## Introduction
A workflow to create and automatically analyze short Molecular Dynamics (MD) simulation from a set of drug and excipient SMILES strings to determine potential co-aggregation.

Dependency
- Python (3.6>)
- Python Packages : Openbabel, Ambertools, MDTraj, Parmed, OpenMM (7.2>)
- Software (packmol)
- TVRegDiff (Github Repository)

## Order of pipeline
0. Prepare 
	- a tab-delimited file containing list of molecules with SMILES per line (e.g. drug_smiles.txt)
	- a tab-delimited file containing component and number of each molecules to simulate (e.g. pair_composition.txt)
1. From SMILES list, create a 3D conformation in Tripos Mol2 format.
2. Create 
	- a simulation box in PDB format
	- a parameter and topology file in prmtop format (AMBER)
	- and a coordinate file in inpcrd format (AMBER)
3. Edit and launch run_OpenMM.sh file with wanted simulation parameters.
	Note: Refer to OpenMM documentation for detailed simulation parameters.
	The code returns a text file called 'simfile', which can be used to check the queue of the simulations.
4. (Optional) Analyze the simulated drug-excipient pairs.
	Note: The code requires mol2, frcmod, and completed simulation file.

### A Complete Example
Step 0. Download sample data.

	
	git clone https://github.com/DanReker/CoAggregators.git
	

Step 1. Create a 3D conformation in Tripos Mol2 formatted files.

	
	./make_mol.sh ../data/selected_drug_smiles.txt
	./make_mol.sh ../data/selected_excipient_smiles.txt
	
	
Step 2. Create all required pre-simulation files.

	
	python ./make_box.py -i ../data/pair_composition.txt
	

Step 3. Run simulations

	
	./run_OpenMM.sh simulation_directory ../data/pair_composition.tsv
	
where simulation_directory links to a directory where simulation files will be stored

Step 4. Automated analysis of simulations 
```
python ./analysis.py -i pair_composition.tsv -i outfile.csv -d simulation_directory
```
