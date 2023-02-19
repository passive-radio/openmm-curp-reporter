# Mockup of Python Reporter class which reports interatomic quantity

mockup python file: [InterAtomicReporter.py](./InterAtomicReporter.py)

test of md simulation with InterAtomicReporter reporter: [test/try_mockup.py](./test/test_interatomic_reporter.py)

Goal image
```python
# Input Files
pdb = PDBFile('input/7r98-processed.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Some Configuration
# ...

# Prepare the Simulation
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
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

# Simulate
simulation.reporters.append(HeatFluxReporter(outpath="out/heat_flux.nc", interval=1000, decomp=True, cutoff=1.0*nanometer))
simulation.currentStep = 0
simulation.step(steps)
```