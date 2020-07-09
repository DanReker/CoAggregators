'''
Generate a simulation pdb box and AMBER-recognized input files from a tab-delimited file that contains the molecule components and the number of each molecules per line.
Example line entry: 'Sorafenib  indomethacin    2   2'

usage: make_box.py [-h] -i INFILE [-s] [-b BUFFER] [-# SEED]
	optional arguments:
	  -h, --help            show this help message and exit
	  -i INFILE, --infile INFILE
				Input file that contains names and number of molecules
	  -s, --solvate         If activated, creates a model solvated in explicit water (TIP3P)
	  -b BUFFER, --buffer BUFFER
				Buffering distance between waterbox and molecule box
				(default=10 Angstroms)
	  -# SEED, --seed SEED  If specified, creates a configuration with given seed
				number (default=random)
'''


import mdtraj as md
import parmed
import pandas as pd
import numpy as np
import os
from mdtraj.utils.delay_import import import_
from distutils.spawn import find_executable
import logging
import subprocess
import sys
import argparse

logger = logging.getLogger(__name__)


def getoutput(cmd):
    """Compatibility function to substitute deprecated commands.getoutput in Python2.7 (Original code : openmoltools.amber)""" 
    try:
        out = subprocess.getoutput(cmd)
    except AttributeError:
        out = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT,
                               stdout=subprocess.PIPE).communicate()[0]
    try:
        return str(out.decode())
    except:
        return str(out)

def mol2_to_pdb(mol2_filename):
	#convert PDB file with this name, to be used for Packmol input
	struct = md.load_mol2(mol2_filename)
	struct.save_pdb(mol2_filename[:-4]+'pdb')

def approximate_volume(pdb_filenames, n_molecules_list, box_scaleup_factor=2.0):
    """Approximate the appropriate box size based on the number and types of atoms present. (Original code : openmoltools.packmol)
    Parameters
    ----------
    pdb_filenames : list(str)
        List of pdb filenames for each component of mixture.
    n_molecules_list : list(int)
        The number of molecules of each mixture component.
    box_scaleup_factor : float, optional, default = 2.0
        Factor by which the estimated box size is increased
    Returns
    -------
    box_size : float
        The size of the box to generate.  In ANGSTROMS.
    Notes
    -----
    By default, boxes are very large for increased stability, and therefore may 
    require extra time for energy minimization and equilibration.
    """
    volume = 0.0  # in cubic angstroms
    for k, (pdb_file) in enumerate(pdb_filenames):
        molecule_volume = 0.0
        molecule_trj = md.load(pdb_filenames[k])
        for atom in molecule_trj.topology.atoms:
            if atom.element.symbol == 'H':
                molecule_volume += 5.0  # approximated from bondi radius = 1.06 angstroms
            else:
                molecule_volume += 15.0  # approximated from bondi radius of carbon = 1.53 angstroms
        volume += molecule_volume * n_molecules_list[k]
    box_size = volume**(1.0/3.0) * box_scaleup_factor
    return box_size

PACKMOL_PATH = find_executable("packmol")


HEADER_TEMPLATE = """
# Mixture 
tolerance %f
filetype pdb
output %s
add_amber_ter
seed %d
"""

BOX_TEMPLATE = """
structure %s
  number %d 
  inside box 0. 0. 0. %f %f %f
end structure
"""

TLEAP_TEMPLATE = """
source leaprc.gaff
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
addions box Na+ 0
addions box Cl- 0
setbox box centers
savepdb box %(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

TLEAP_TEMPLATE_SOLV = """
source leaprc.gaff
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
addions box Na+ 0
addions box Cl- 0
solvatebox box TIP3PBOX %(buffer_distance)d
setbox box centers
savepdb box solv_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""


def pack_box(pdb_filenames, n_molecules_list, tolerance=2.0, box_size=None, seed=-1):
    """Run packmol to generate a box containing a mixture of molecules. (Original code : openmoltools.packmol)
    Parameters
    ----------
    pdb_filenames_or_trajectories : list({str, Trajectory})
        List of pdb filenames or trajectories for each component of mixture.
        If this is a list of trajectories, the trajectories will be saved to
        as temporary files to be run in packmol. Water molecules must have
        MDTraj-standard residue and atom names as defined in
        mdtraj/formats/pdb/data/pdbNames.xml.
    n_molecules_list : list(int)
        The number of molecules of each mixture component.
    tolerance : float, optional, default=2.0
        The minimum spacing between molecules during packing.  In ANGSTROMS!
    box_size : float, optional
        The size of the box to generate.  In ANGSTROMS.
        Default generates boxes that are very large for increased stability.
        May require extra time for energy minimization and equilibration.
    seed : default = -1 (random). User may specify seed number for consistent configuration generation
    Returns
    -------
    trj : MDTraj.Trajectory
        Single frame trajectory with mixture box.
    Notes
    -----
    Water molecules must have MDTraj-standard residue and atom names as defined
    in mdtraj/formats/pdb/data/pdbNames.xml, otherwise MDTraj won't be able to
    perceive the bonds and the Topology of the returned Trajectory will be incorrect.
    Be aware that MDTraj uses nanometers internally, but packmol uses angstrom
    units. The present function takes `tolerance` and `box_size` in angstrom
    units, but the output trajectory will have data in nm.
    Also note that OpenMM is pretty picky about the format of unit cell input, 
    so use the example in tests/test_packmol.py to ensure that you do the right thing.
    See Also
    --------
    standardize_water
        Standardize residue and atom names of a water molecule.
    """
    assert len(pdb_filenames) == len(
    	n_molecules_list), "Must input same number of pdb filenames as num molecules"

    if PACKMOL_PATH is None:
        raise(IOError("Packmol not found, cannot run pack_box()"))

    trj_i=[]
    for obj in pdb_filenames:
        trj_i.append(md.load(obj))
            
    # Approximating volume to initialize box.
    if box_size is None:
        box_size = approximate_volume(pdb_filenames, n_molecules_list)

	# Adjust box_size for periodic box. Packmol does not explicitly
	# support periodic boundary conditions and the suggestion on
	# their docs is to pack in a box 2 angstroms smaller. See
	# http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml#pbc
    packmol_box_size = box_size - 2  # angstroms

	# The path to packmol's output PDB file. Concatenate name of molecules
    output_filename = ''
    for i in range(len(pdb_filenames)-1):	
        output_filename = output_filename + pdb_filenames[i][:-4]+'_'
    output_filename = 'box_'+ output_filename + pdb_filenames[len(pdb_filenames)-1][:-4]+'.pdb' 

	# Create input file for packmol.
    header = HEADER_TEMPLATE % (tolerance, output_filename, seed)
    for k in range(len(pdb_filenames)):
        filename = pdb_filenames[k]
        n_molecules = n_molecules_list[k]
        header += BOX_TEMPLATE % (filename, n_molecules, packmol_box_size, packmol_box_size, packmol_box_size)
    print(header)
    packmol_filename = output_filename+'.inp'
    with open(packmol_filename, 'w') as file_handle:
        file_handle.write(header)

	# Run packmol and load output PDB file.
    os.system("%s < %s" % (PACKMOL_PATH, file_handle.name))
    trj = md.load(output_filename)

    assert trj.topology.n_chains == sum(n_molecules_list), "Packmol error: molecules missing from output"

    #Begin hack to introduce bonds for the MISSING CONNECT ENTRIES THAT PACKMOL FAILS TO WRITE

    top, bonds = trj.top.to_dataframe()
    bonds_i = [t.top.to_dataframe()[1] for t in trj_i]

    offset = 0
    bonds = []
    for i in range(len(trj_i)):
        n_atoms = trj_i[i].n_atoms
        for j in range(n_molecules_list[i]):
            list(bonds).extend(bonds_i[i] + offset)
            offset += n_atoms

    bonds = np.array(bonds)
    trj.top = md.Topology.from_dataframe(top, bonds)

    # Set the requested box size.
    trj.unitcell_vectors = np.array([np.eye(3)]) * box_size / 10.

    return output_filename


def build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename, solvation=False, buffer_distance=10):
    """Create a prmtop and inpcrd from a collection of mol2 and frcmod files
    as well as a single box PDB.  We have used this for setting up
    simulations of binary mixtures. (Original code : openmoltools.amber)
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
    solvation : Boolean, optional. Default: False
        Boolean for whether the system should be solvated explicitly or not. If true, the system will be solvated in TIP3P water model using tleap
    buffer_distance : int, optional. Default: 10
        If solvation is true, will add water molecules with a buffering distance of 10 Angstrom unless specified otherwise.
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

    #Build absolute paths of output files so we can copy them back
    prmtop_filename = os.path.abspath(prmtop_filename)
    inpcrd_filename = os.path.abspath(inpcrd_filename)

    all_names = [md.load(filename).top.residue(0).name for filename in mol2_filenames]

    mol2_section = "\n".join("%s = loadmol2 %s" % (
        all_names[k], filename) for k, filename in enumerate(mol2_filenames))
    #If non-GAFF water is present, load desired parameters for that water as well.
    amberparams_section = "\n".join("loadamberparams %s" % (filename) for k, filename in enumerate(frcmod_filenames))

    if solvation == True:
        tleap_commands = TLEAP_TEMPLATE_SOLV % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, buffer_distance=buffer_distance,
                                            box_filename=box_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)
    elif solvation == False:
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


# MAIN CODE	
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-i","--infile", type=open, help="Input file that contains names and number of molecules", required=True)
parser.add_argument("-s", "--solvate", action="store_true", help="If activated, creates a model solvated in explicit water (TIP3P)")
parser.add_argument("-b", "--buffer", type=int, default=10, help="Buffering distance between waterbox and molecule box (default=10 Angstroms)")
parser.add_argument("-#", "--seed", type=int, default=-1, help="If specified, creates a configuration with given seed number (default=random)")
args = parser.parse_args()

pair_df = pd.read_csv(args.infile, sep='\t',header=None)
splitpoint = len(pair_df.columns)/2
mol_df = pair_df.loc[:,0:splitpoint-1]
num_df = pair_df.loc[:,splitpoint:]

if args.solvate:
    solvation = True
    buffer_distance = args.buffer
else:
    solvation = False
    buffer_distance = None

seed=args.seed

for i,row in pair_df.iterrows():
    pdb_list = []
    num_list = []
    mol2_filenames = []
    frcmod_filenames = []
    mix_name = ''
    for mol in mol_df.loc[i]:
        mol2_to_pdb(mol+'.mol2')
        pdb_list += [mol+'.pdb']
        mol2_filenames += [mol+'.mol2']
        frcmod_filenames += [mol+'.frcmod']
        mix_name = mix_name + mol + '_'
    for num in num_df.loc[i]:
        num_list += [int(num)]

    final_name = mix_name[:-1]
    box_filename = pack_box(pdb_list, num_list, seed=seed)
    tleap_cmd = build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, final_name+'.prmtop', final_name+'.inpcrd', solvation=solvation, buffer_distance=buffer_distance)

    file_handle = open(final_name+".leap.in", 'w')
    file_handle.writelines(tleap_cmd)
    file_handle.close()