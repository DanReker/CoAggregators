import os
from distutils.spawn import find_executable
import mdtraj as md
import numpy as np
import pandas as pd
import parmed
from tvregdiff import TVRegDiff
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmoltools.utils import getoutput
import argparse
import logging

'''
Generates a dataframe containing parameters for drug-excipient pair used for machine learning post-analysis. The full dataframe includes 19 different values for each simulation pair

usage: python analysis.py [-h] -i INFILE -o OUTFILE -d DIR
example: python ./analysis.py -i ../data/pair_composition.tsv -o outfile.csv -d simulation_directory
	optional arguments:
	  -h, --help            show this help message and exit
	  -i INFILE, --infile INFILE
				Input file that contains names and number of molecules
	  -o OUTFILE, --outfile OUTFILE
				Name of the csv file
	  -d DIR, --directory DIR
				Provide the directory in which simulation files are located             
'''

def build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename):
    """Create a prmtop and inpcrd from a collection of mol2 and frcmod files as well as a single box PDB.  
    We have used this for setting up simulations of binary mixtures.
    - Original code by : Chodera Lab / Openmoltools project (https://github.com/choderalab/openmoltools)

    Parameters
    ----------
    mol2_filenames : list(str)
        Filenames of GAFF flavored mol2 files.  Each must contain exactly
        ONE ligand.
    frcmod_filenames : str
        Filename of input GAFF frcmod filenames.
    box_filename : str
        Filename of PDB containing an arbitrary box of the mol2 molecules.
    prmtop_filename : str
        output prmtop filename.  Should have suffix .prmtop
    inpcrd_filename : str
        output inpcrd filename.  Should have suffix .inpcrd
    water_model : str, optional. Default: "TIP3P"
        String specifying water model to be used IF water is present as a component of the mixture. Valid options are currently "TIP3P", "SPC", or None. If None is specified, flexible GAFF-water will be used as for any other solute (old behavior).
    Returns
    -------
    tleap_commands : str
        The string of commands piped to tleap for building the prmtop
        and inpcrd files.  This will *already* have been run, but the
        output can be useful for debugging or archival purposes. However,
        this will reflect temporary file names for both input and output
        file as these are used to avoid tleap filename restrictions.
    Notes
    -----
    This can be easily broken if there are missing, duplicated, or
    inconsistent ligand residue names in your box, mol2, and frcmod files.
    You can use mdtraj to edit the residue names with something like
    this: trj.top.residue(0).name = "L1"
    """

    # Check for one residue name per mol2 file and uniqueness between all mol2 files
    all_names = set()
    for filename in mol2_filenames:
        t = md.load(filename)
        names = set([r.name for r in t.top.residues])

        if len(names) != 1:
            raise(ValueError("Must have a SINGLE residue name in each mol2 file."))

        all_names = all_names.union(list(names))

    if len(all_names) != len(mol2_filenames):
        raise(ValueError("Must have UNIQUE residue names in each mol2 file."))
    if len(mol2_filenames) != len(frcmod_filenames):
        raise(ValueError("Must provide an equal number of frcmod and mol2 file names."))

    #Get number of files
    nfiles = len(mol2_filenames)

    #Build absolute paths of input files so we can use context and temporary directory
    infiles = mol2_filenames + frcmod_filenames + [box_filename]
    infiles = [os.path.abspath(filenm) for filenm in infiles]

    #Build absolute paths of output files so we can copy them back
    prmtop_filename = os.path.abspath(prmtop_filename)
    inpcrd_filename = os.path.abspath(inpcrd_filename)

    #Use temporary directory and do the setup
    with md.utils.enter_temp_directory():
        all_names = [md.load(filename).top.residue(
            0).name for filename in mol2_filenames]

        mol2_section = "\n".join("%s = loadmol2 %s" % (
            all_names[k], filename) for k, filename in enumerate(mol2_filenames))

        amberparams_section = "\n".join("loadamberparams %s" % (filename) for k, filename in enumerate(frcmod_filenames))

        tleap_commands = TLEAP_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section,
                                               box_filename=box_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)
        print(tleap_commands)

        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleap_commands)
        file_handle.close()

        logger.debug('Running tleap in temporary directory.')
        cmd = "tleap -f %s " % file_handle.name
        logger.debug(cmd)

        output = getoutput(cmd)
        logger.debug(output)

    return tleap_commands

def find_gaff_dat():
    AMBERHOME = None

    try:
        AMBERHOME = os.environ['AMBERHOME']
    except KeyError:
        pass

    if AMBERHOME is None:
        full_path = find_executable("parmchk2")
        try:
            AMBERHOME = os.path.split(full_path)[0]
            AMBERHOME = os.path.join(AMBERHOME, "../")
        except:
            raise(ValueError("Cannot find AMBER GAFF"))

    if AMBERHOME is None:
        raise(ValueError("Cannot find AMBER GAFF"))

    return os.path.join(AMBERHOME, 'dat', 'leap', 'parm', 'gaff.dat')

#########################################################################################################

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-i","--infile", type=open, help="Input file that contains name(s) of simulation file(s).", required=True)
parser.add_argument("-o","--outfile", type=str, default='outfile.csv', help="Specify the name of the output csv file.", required=True)
parser.add_argument("-d","--directory", type=str, help="Provide the directory in which simulation files are located.", required=True)
args = parser.parse_args()

# simulation file names
simfiles = []
drugs = []
excip = []
for line in args.infile:
    simfiles += ['openmm_'+line.split("\t")[0].strip()+'_'+line.split("\t")[1].strip()+'.pdb']
    drugs += [line.split("\t")[0].strip()]
    excip += [line.split("\t")[1].strip()]


names = pd.DataFrame(data=None, index=None, columns=["sim_filename", "drug_name", "excip_name", "d_mol2_filename", "d_frcmod_filename","e_mol2_filename", "e_frcmod_filename"])
names['sim_filename'] = simfiles
names['drug_name'] = drugs
names["d_mol2_filename"] = names["drug_name"] + ".mol2"
names["d_frcmod_filename"] = names["drug_name"] + ".frcmod"
names['excip_name'] = excip
names["e_mol2_filename"] = names["excip_name"] + ".mol2"
names["e_frcmod_filename"] = names["excip_name"] + ".frcmod"

# Start with a dataframe that contains column of simulation filenames
count = 0
# 1~4 was selected by order of drug-excipient residue pair with smallest average distance to largest
df = pd.DataFrame(data=None, index=None, columns=["drug", "excipient", "contact_av", "LF_KE", "LF_PE", "AvDist_1", "AvDist_2", "AvDist_3", "AvDist_4", "d(distance)_max_1", "d(distance)_max_2", "d(distance)_max_3", "d(distance)_max_4",  "d(distance)_avg_1", "d(distance)_avg_2", "d(distance)_avg_3", "d(distance)_avg_4", "Variance_1", "Variance_2", "Variance_3", "Variance_4"])
df['drug'] = drugs
df['excipient'] = excip

logger = logging.getLogger(__name__)

TLEAP_TEMPLATE = """
source leaprc.gaff
source oldff/leaprc.ff99SB
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
setbox box centers
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

PAIR_DIR=os.path.abspath(args.directory)+'/'

for i in df.index: # distance_df.index:
    simfile = names.loc[i, "sim_filename"]
    try:
        print(simfile+" loading...\n")
        traj = md.load_pdb(PAIR_DIR+simfile)
        print(simfile+" loaded\n")
        # shorten simulation files to 2000 frames
        if traj.n_frames == 10000:
            traj = traj[::5]
        if traj.n_frames > 10000:
            scale = traj.n_frames/2000
            traj = traj[::scale]

        # Create last frame snapshot for later analysis
        last_frame = traj.slice(traj.n_frames-1)
        lf_file = 'lastFrame_' + simfile
        last_frame.save_pdb(PAIR_DIR+lf_file)

        drug_molecule_filename = PAIR_DIR+ names.loc[i,"d_mol2_filename"]
        drug_frcmod_filename = PAIR_DIR + names.loc[i, "d_frcmod_filename"]
        excipient_molecule_filename = PAIR_DIR + names.loc[i,"e_mol2_filename"]
        excipient_frcmod_filename =PAIR_DIR + names.loc[i,"e_frcmod_filename"]

        prmtop_filename = PAIR_DIR + 'LF_'+ names.loc[i,"drug_name"]+"_"+names.loc[i,"excip_name"]+".prmtop"
        inpcrd_filename = PAIR_DIR + 'LF_'+ names.loc[i,"drug_name"]+"_"+names.loc[i,"excip_name"]+".inpcrd"

        # reset the residue name according to lastframe file
        box_filename = PAIR_DIR + 'lastFrame_openmm_' + names.loc[i, "drug_name"]+"_"+names.loc[i, "excip_name"]+".pdb"
        lf_box = md.load_pdb(box_filename)
        drg_res = str(lf_box.topology.residue(0))[:3]
        exc_res = str(lf_box.topology.residue(2))[:3]

        drg_struct = parmed.load_file(drug_molecule_filename)
        drg_struct.name = drg_res
        drg_mol2file = parmed.formats.Mol2File
        drg_mol2file.write(drg_struct, drug_molecule_filename)

        exc_struct = parmed.load_file(excipient_molecule_filename)
        exc_struct.name = exc_res
        exc_mol2file = parmed.formats.Mol2File
        exc_mol2file.write(exc_struct, excipient_molecule_filename)

        mol2_filenames = [drug_molecule_filename, excipient_molecule_filename]
        frcmod_filenames = [drug_frcmod_filename, excipient_frcmod_filename]

        #Calculation starts
        distance = md.compute_contacts(traj, contacts=[[0,2],[0,3],[1,2],[1,3]], scheme='closest-heavy', periodic=False)

        # contact number (two residues closer than 4.5 Angstrom)
        total_contact = 0
        average_contact = 0
        for j in range(len(distance[0])):
            total_contact += sum(distance[0][j,:] < 0.45)
            average_contact = total_contact/len(distance[0])
        df.loc[i, "contact_av"] = average_contact

        # get average distance and order from least to greatest for all residue pairs
        AveDist = np.sum(distance[0][:, :], axis=0)/len(distance[0])
        dist_order = np.argsort(AveDist)

        ### Average distance
        df.loc[i, "AvDist_1"] = AveDist[dist_order[0]]
        df.loc[i, "AvDist_2"] = AveDist[dist_order[1]]
        df.loc[i, "AvDist_3"] = AveDist[dist_order[2]]
        df.loc[i, "AvDist_4"] = AveDist[dist_order[3]]

        # subset trajectory for the pair with minAveDist
        traj_pair_1 = distance[0][:, dist_order[0]]
        traj_pair_2 = distance[0][:, dist_order[1]]
        traj_pair_3 = distance[0][:, dist_order[2]]
        traj_pair_4 = distance[0][:, dist_order[3]]

        print("Calculating derivative...")
        d1 = TVRegDiff(traj_pair_1, 10, 1, scale='large', plotflag=0, diagflag=0)
        d2 = TVRegDiff(traj_pair_2, 10, 1, scale='large', plotflag=0, diagflag=0)
        d3 = TVRegDiff(traj_pair_3, 10, 1, scale='large', plotflag=0, diagflag=0)
        d4 = TVRegDiff(traj_pair_4, 10, 1, scale='large', plotflag=0, diagflag=0)

        df.loc[i, "d(distance)_max_1"] = d1.max()
        df.loc[i, "d(distance)_max_2"] = d2.max()
        df.loc[i, "d(distance)_max_3"] = d3.max()
        df.loc[i, "d(distance)_max_4"] = d4.max()

        df.loc[i, "d(distance)_avg_1"] = d1.mean()
        df.loc[i, "d(distance)_avg_2"] = d2.mean()
        df.loc[i, "d(distance)_avg_3"] = d3.mean()
        df.loc[i, "d(distance)_avg_4"] = d4.mean()

        df.loc[i, "Variance_1"] = np.var(traj_pair_1)
        df.loc[i, "Variance_2"] = np.var(traj_pair_2)
        df.loc[i, "Variance_3"] = np.var(traj_pair_3)
        df.loc[i, "Variance_4"] = np.var(traj_pair_4)
        
        tleap_cmd = build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename)
        tleap_filename = PAIR_DIR + names.loc[i, "drug_name"]+"_"+names.loc[i, "excip_name"]+".leap.in"
        file_handle = open(tleap_filename, 'w')
        file_handle.writelines(tleap_cmd)
        file_handle.close()
            
        print("Calculating energy...")
        # load last frame and setup simulation and integrator. get the potential/kinetic energy state at the last frame
        prmtop = AmberPrmtopFile(PAIR_DIR+'LF_'+simfile[7:-4]+'.prmtop')
        inpcrd = AmberInpcrdFile(PAIR_DIR+'LF_'+simfile[7:-4]+'.inpcrd')
        system = prmtop.createSystem(implicitSolvent=OBC2, nonbondedMethod=NoCutoff, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(prmtop.topology, system, integrator)
        simulation.context.setPositions(inpcrd.positions)
        PE = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        df.loc[i, 'LF_PE'] = PE.value_in_unit(unit=kilojoule/mole)
        KE = simulation.context.getState(getEnergy=True).getKineticEnergy()
        df.loc[i, 'LF_KE'] = KE.value_in_unit(unit=kilojoule/mole)

        count = count + 1
        print("Completed "+ str(count) + "\n")      
    except Exception as e:
        print(simfile+" skipped!\n")
        print(e)
        continue

df.to_csv(PAIR_DIR+args.outfile, sep='\t', index=False)
