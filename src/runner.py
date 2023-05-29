#!/usr/bin/env python
from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm
from sys import stdout
import sys
from myreporter import ForceReporter
import argparse

def save_lastsnapshot(simulation, outname='output'):
    print('Saving the last snapshot...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f"{outname}.pdb", 'w'))

class Runner():
    def __init__(self, filename) -> None:
        self.check_openmm_version()

        print('Loading...')
        self.pdb = PDBFile(filename)
        self.forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        self.modeller = Modeller(self.pdb.topology, self.pdb.positions)
        self.watermodel = 'tip3p'

    @staticmethod
    def check_openmm_version():
        required_version = "8.0"
        if openmm.__version__ != required_version:
            sys.exit(f"You are using a version the author did not use.: {openmm.__version__}")

    def prepare_system(self):    
        print('Adding hydrogens...')
        self.modeller.addHydrogens(self.forcefield, pH=7.0)
        
        print('Adding solvent...')
        self.modeller.addSolvent(forcefield=self.forcefield, 
                                 model=self.watermodel, 
                                 padding=1*nanometer, 
                                 ionicStrength=1*molar,
                                 neutralize=True)
                                    #boxShape='dodecahedron' # over openmm 8.0 
        
        print('System building...')
        self.system = self.forcefield.createSystem(self.modeller.topology,
                                                   nonbondedMethod=PME, 
                                                   nonbondedCutoff=1*nanometer,
                                                   constraints=HBonds)
        
        self.integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        # ^ integrator must be defined before creating simulation object that takes it.
        #   For clarity, I redefine an integrator in the function of `equilibriate` because thermostat is used at this point.
        self.simulation = Simulation(self.modeller.topology, self.system, self.integrator)
        self.simulation.context.setPositions(self.modeller.positions)  

    def add_reporters(self, mdsteps=1000, logperiod=10, dcdperiod=10, forceperiod=10):
        print("Setting reporters...")
        ## For standard output
        self.simulation.reporters.append(StateDataReporter(stdout, 
                                                    reportInterval=10*logperiod, 
                                                    step=True,
                                                    time=True, 
                                                    progress=True, 
                                                    remainingTime=True, 
                                                    speed=True, 
                                                    totalSteps=mdsteps, 
                                                    separator='\t')
                                    )
        ## For writing a simulation log file.
        # TODO: add pressure output
        self.simulation.reporters.append(StateDataReporter("md.csv", 
                                                    reportInterval=logperiod, 
                                                    time=True,
                                                    totalEnergy=True, 
                                                    kineticEnergy=True, 
                                                    potentialEnergy=True, 
                                                    temperature=True, 
                                                    density=True, 
                                                    volume=True,
                                                    separator=","
                                                    )
                                    )
        self.simulation.reporters.append(DCDReporter('traj.dcd', dcdperiod))
        self.simulation.reporters.append(CheckpointReporter('checkpoint.chk', logperiod, writeState=False)) 
        # ^ .chk is binary, so specific to the envirinment where it is created. 
        #self.simulation.reporters.append(ForceReporter('force.dat', forceperiod))

    def minimise(self):
        print('Minimizing...')
        self.simulation.minimizeEnergy(maxIterations=100) #can't I report data for em?
        minpositions = self.simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(self.simulation.topology, minpositions, open('min.pdb', 'w'))

    def equilibriate(self):    
        print('Equilibration step via NVT...')
        nvtmdsteps=100
        self.integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        self.simulation.step(nvtmdsteps)
        save_lastsnapshot(self.simulation, "nvt_eq")
        
        # TODO: I should check the volume and pressure to confirm NPT actually works. 
        print('Equilibration step via NPT...')
        nptmdsteps=100
        self.barostat = MonteCarloBarostat(1.0*bar, 300.0*kelvin, 25) 
        self.system.addForce(self.barostat)
        self.simulation.loadCheckpoint('checkpoint.chk')
        self.simulation.step(nptmdsteps)
        save_lastsnapshot(self.simulation, "npt_eq")

    #TODO: abstract this function for NPT NVT ensemble. 
    def nvt_production(self, nsteps=1000, is_xml_output=True):
        print('Production run...')
        if is_xml_output: # We often do not want all the data that is too large, so this option turns off the full output.
            self.simulation.saveState('output.xml')
            
        #self.simulation.loadCheckpoint('checkpoint.chk')
        self.simulation.step(nsteps)  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run molecular dynamics simulations using OpenMM.')
    parser.add_argument('filename', type=str, help='Input PDB file')
    args = parser.parse_args()
    
    runner = Runner(args.filename)
    runner.prepare_system()
    runner.add_reporters()
    runner.minimise()
    runner.equilibriate()
    runner.nvt_production()
