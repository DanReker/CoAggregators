'''
usage: OpenMM.py [-h] -top TOPOLOGY -crd COORDINATE -o OUTPUT [--savechk]
         [--loadchk] [-c {0,1,2,3}] [-p] [-i]
         [--solvent_type {1,2,3,4,5}] [-T TEMPERATURE]
         [--timestep TIMESTEP] [--interval INTERVAL]
         [-l SIMULATION_LENGTH] [--timeit]

    optional arguments:
      -h, --help            show this help message and exit
      -top TOPOLOGY, --topology TOPOLOGY
                Input prmtop file
      -crd COORDINATE, --coordinate COORDINATE
                Input inpcrd file
      -o OUTPUT, --output OUTPUT
                Output pdb filename
      --savechk             If activated, save a checkpoint state at every 10
                intervals
      --loadchk             If activated, load a checkpoint state with .chk
                extension of the output filename
      -c {0,1,2,3}, --constraint {0,1,2,3}
                0 = None (default); 1 = HBonds ; 2 = AllBonds ; 3 =
                HAngles
      -p, --periodic        If activated, runs the simulation in periodic box with
                PME method used for long range interaction (default =
                NoCutoff)
      -i, --implicit        If activated, runs the simulation in implicit water
                solvent (default = vacuum, unless explicitly solvated)
      --solvent_type {1,2,3,4,5}
                1 = HCT (default); 2 = OBC1 ; 3 = OBC2 ; 4 = GBn ; 5 =
                GBn2
      -T TEMPERATURE, --temperature TEMPERATURE
                Set simulation temperature (default = 300K)
      --timestep TIMESTEP   Set simulation time step in units of picosecond
                (default = 0.002 picosecond)
      --interval INTERVAL   Set interval of saved frames in the unit of picosecond
                (default = 10 ps)
      -l SIMULATION_LENGTH, --simulation_length SIMULATION_LENGTH
                Set duration of simulation time in units of nanosecond
                (default = 20ns)
      --timeit              If activated, creates time_log file
'''

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import time
import os
import argparse

CUR_DIR = os.getcwd()
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("-top","--topology", help="Input prmtop file", required=True)
parser.add_argument("-crd","--coordinate", help="Input inpcrd file", required=True)
parser.add_argument("-o","--output", help="Output pdb filename", required=True)
parser.add_argument("--savechk", action="store_true", help="If activated, save a checkpoint state at every 10 intervals")
parser.add_argument("--loadchk", action="store_true", help="If activated, load a checkpoint state with .chk extension of the output filename")
parser.add_argument('-c','--constraint', type=int, default=0, choices=range(0, 4), help="0 = None (default); 1 = HBonds ; 2 = AllBonds ; 3 = HAngles")
parser.add_argument('-p','--periodic', action="store_true", help="If activated, runs the simulation in periodic box with PME method used for long range interaction (default = NoCutoff)")
parser.add_argument('-i','--implicit', action="store_true", help="If activated, runs the simulation in implicit water solvent (default = vacuum, unless explicitly solvated)")
parser.add_argument('--solvent_type', type=int, default=1, choices=range(1, 6), help="Actiavted if '-i' is specified. 1 = HCT ; 2 = OBC1 ; 3 = OBC2 (default) ; 4 = GBn ; 5 = GBn2 ")
parser.add_argument('-T','--temperature', type=int, default=300, help="Set simulation temperature (default = 300K)")
parser.add_argument('--timestep', type=float, default=0.002, help="Set simulation time step in units of picosecond (default = 0.002 picosecond)")
parser.add_argument('--interval', type=int, default=10, help="Set interval of saved frames in the unit of picosecond (default = 10 ps)")
parser.add_argument('-l','--simulation_length', type=int, default=20, help="Set duration of simulation time in units of nanosecond (default = 20ns)")
parser.add_argument('--timeit', action="store_true", help="If activated, creates time_log file")
args = parser.parse_args()

constraint_level=None
if args.constraint==1:
    constraint_level="HBonds"
elif args.constraint==2:
    constraint_level="AllBonds"
elif args.constraint==3:
    constraint_level="HAngles"

if args.periodic:
    periodicity = PME
else:
    periodicity = NoCutoff

if not args.implicit:
    solvent_type=None #0
else:
    if args.solvent_type==1:
        solvent_type=HCT
    if args.solvent_type==2:
        solvent_type=OBC1
    if args.solvent_type==3:
        solvent_type=OBC2
    if args.solvent_type==4:
        solvent_type=GBn
    if args.solvent_type==5:
        solvent_type=GBn2        

interval = int(args.interval / args.timestep)
simlength = int(args.simulation_length * 1000 / args.timestep)


### Actual code
try:
    prmtop = AmberPrmtopFile(args.topology)
    inpcrd = AmberInpcrdFile(args.coordinate)
    system = prmtop.createSystem(nonbondedMethod=periodicity, constraints=constraint_level, implicitSolvent=solvent_type)
    integrator = LangevinIntegrator(args.temperature*kelvin, 1/picosecond, args.timestep*picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    simulation.minimizeEnergy()
    if args.loadchk:
        simulation.loadCheckpoint(args.output+'.chk')
    simulation.reporters.append(PDBReporter(args.output, interval))
    simulation.reporters.append(StateDataReporter(args.output+'.log', interval, step=True,
            potentialEnergy=True, temperature=True))
    if args.savechk:
        simulation.reporters.append(CheckpointReporter(args.output+'.chk', 10*interval))
    simulation.step(simlength)
    # create time_log files
    if args.timeit:
        log_filename = CUR_DIR + args.output + ".time_log"
        file_handle = open(log_filename, 'w')
        file_handle.writelines("--- %s seconds ---" % (time.time() - start_time))
        file_handle.close()
except:
    failfile = open("failed.txt", "a")
    failfile.writelines(args.output + '\n')
    failfile.close()
