import unittest
from openmm.app import *
from openmm import Platform, MonteCarloBarostat, LangevinMiddleIntegrator 
import openmm as mm
from openmm.unit import *

import sys
sys.path.append('../')
from reporter import InterAtomicReporter

class ReporterTest(unittest.TestCase):
    def setUp(self) -> None:
        return super().setUp()
    
    def tearDown(self) -> None:
        return super().tearDown()

    def test_reporter(self):
        # Input Files
        pdb = PDBFile('../input/7r98-processed.pdb')
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

        # System Configuration
        nonbondedMethod = PME
        nonbondedCutoff = 1.0*nanometers
        ewaldErrorTolerance = 0.0005
        constraints = HBonds
        rigidWater = True
        constraintTolerance = 0.000001

        # Integration Options
        dt = 0.002*picoseconds
        temperature = 300*kelvin
        friction = 1.0/picosecond
        pressure = 1.0*atmospheres
        barostatInterval = 25

        # Simulation Options
        steps = 10000
        equilibrationSteps = 5000
        platform = Platform.getPlatformByName('CUDA')
        platformProperties = {'Precision': 'single'}

        # Prepare the Simulation
        print('Building system...')
        topology = pdb.topology
        positions = pdb.positions
        system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
            constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
        system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
        integrator = LangevinMiddleIntegrator(temperature, friction, dt)
        integrator.setConstraintTolerance(constraintTolerance)
        simulation = Simulation(topology, system, integrator, platform, platformProperties)
        simulation.context.setPositions(positions)

        # Minimize and Equilibrate
        print('Performing energy minimization...')
        simulation.minimizeEnergy()
        print('Equilibrating...')
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(equilibrationSteps)

        # Simulate
        print('Simulating...')
        simulation.reporters.append(InterAtomicReporter.InterAtomicReporter('log.txt', 100))
        simulation.currentStep = 0
        simulation.step(steps)

if __name__ == '__main__':
    unittest.main()