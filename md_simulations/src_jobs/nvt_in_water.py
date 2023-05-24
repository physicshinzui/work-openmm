#!/usr/bin/env python
from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm
from sys import stdout
import sys

print("openmm version: ", openmm.__version__)
if openmm.__version__ != "8.0":
    sys.exit(f"You are using a version that the author did not use.: {openmm.__version__}")

filename=sys.argv[1]
watermodel = 'tip3p'

print('Loading...')
pdb = PDBFile(filename)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)

print('Adding hydrogens...')
modeller.addHydrogens(forcefield, pH=7.0)

print('Adding solvent...')
modeller.addSolvent(forcefield=forcefield, 
                    model=watermodel, 
                    padding=1*nanometer, 
                    ionicStrength=1*molar,
                    neutralize=True)
                    #boxShape='dodecahedron' # over openmm 8.0 

print('System building...')
system = forcefield.createSystem(modeller.topology,
                                nonbondedMethod=PME, 
                                nonbondedCutoff=1*nanometer,
                                constraints=HBonds)

#integrator = VerletIntegrator(0.002*picoseconds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

mdsteps = 1000
logperiod = 10
dcdperiod = 10
print("Setting reporters...")
simulation.reporters.append(StateDataReporter(stdout, 10*logperiod, step=True,
                            time=True, progress=True, remainingTime=True, speed=True, totalSteps=mdsteps, separator='\t'))
simulation.reporters.append(StateDataReporter("md.csv", reportInterval=logperiod, time=True,
                            totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, density=True, volume=True,
                            separator=","))
simulation.reporters.append(DCDReporter('traj.dcd', dcdperiod))
simulation.reporters.append(CheckpointReporter('checkpoint.chk', logperiod, writeState=False)) #.chk is binary, so specific to the envirinment 

print('Minimizing...')
simulation.minimizeEnergy(maxIterations=100) #can't I report data for em?
minpositions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, minpositions, open('min.pdb', 'w'))

nvtmdsteps=1000
print('Equilibration step via NVT...')
simulation.step(nvtmdsteps)

nptmdsteps=500
print('Equilibration step via NPT...')
barostat = MonteCarloBarostat(1.0*bar, 300.0*kelvin, 25) 
system.addForce(barostat)
#simulation.loadCheckpoint('checkpoint.chk')
#simulation.step(nptmdsteps)

print('Saving...')
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))
print('Done')

#print('Restarting...')
#simulation.loadCheckpoint('checkpoint.chk')
#simulation.step(1000)
#system = forcefield.createSystem(modeller.topology, 
#                                 nonbondedMethod=PME, 
#                                 nonbondedCutoff=1*nanometer,
#                                 constraints=HBonds)
#integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
#simulation = Simulation(modeller.topology, system, integrator)
#simulation.context.setPositions(modeller.positions)
#simulation.minimizeEnergy()

#simulation.reporters.append(PDBReporter('output.pdb', 1000))
#simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
#        potentialEnergy=True, temperature=True))
#simulation.step(2000)
