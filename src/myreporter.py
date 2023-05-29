from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm

class ForceReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, True, False, None)

    def report(self, simulation, state):
        forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
        time = state.getTime().value_in_unit(picosecond)
        self._out.write(f"# {time:g} ps \n")
        for f in forces:
            self._out.write(f"{f[0]:g} {f[1]:g} {f[2]:g}\n")

# TODO: 
# Fix this error: openmm.OpenMMException: Invoked getVelocities() on a State which does not contain velocities.
# So at the moment, this class does not handle velocities. 
class XVFReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, True, False, None)

    def report(self, simulation, state):
        positions = state.getPositions().value_in_unit(nanometer)
        velocities = state.getVelocities().value_in_unit(nanometer/picosecond)
        forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
        time = state.getTime().value_in_unit(picosecond)
        self._out.write(f"# {time} ps\n")
        for r, v, f in zip(positions, velocities, forces):
            self._out.write("{r[0]} {r[1]} {r[2]} {v[0]} {v[1]} {v[2]} {f[0]} {f[1]} {f[2]}")