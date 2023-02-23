from __future__ import absolute_import
__author__ = "Yudai Okubo"
__version__ = "0.0.1"

import openmm as mm
from openmm.app import Simulation, Topology
import openmm.unit

import sys

class InterAtomicReporter(object):
    """The base class of Reporter

    Parameters
    ----------
    object : 
        
    """
    def __init__(self, file, reportInterval):
        self.__reportInterval = reportInterval
        
        try:
            self._out = open(file, 'wb')
        except: 
            self._out = open(file,'r+b')
            
        self._out.write('step,donor,accepter,length,quantity\n'.encode('UTF-8'))
        
    def describeNextReport(self, simulation: Simulation):
        steps = self.__reportInterval - simulation.currentStep%self.__reportInterval
        return (steps, self, True, True, True, True)
    
    def report(self, simulation: Simulation, state: mm.State):
        
        """
        # ToDo
        Obtain atoms inside simulation system.
        Get interAtomic quantities from loaded C++ module which runs any interatomic calculation.
        
        ```python
        donor, accepter, quantity = any, any, any
        self._out.write(f'{donor},{accepter},{quantity}')
        ```
        """
        
        system = simulation.system
        simulation_time = state.getTime()
        positions = state.getPositions()
        forces = system.getForces()
        vels = state.getVelocities()
        # params = list(state.getParameters())
        topology = simulation.topology
        bonds = topology.bonds()
        
        for bond in bonds:
            print(bond)
            break
        
        class AtomInfo():
            def __init__(self) -> None:
                self.name = ""
                self.element = None
                self.index = None
                self.residue = None
                self.id = None
            def addNameAndId(self, name, element, index, residue, id):
                self.name = name
                self.element = element
                self.index = index
                self.residue = residue
                self.id = id
            def addNonbondedParams(self, charge, sigma, epsilon):
                self.charge = charge
                self.sigma = sigma
                self.epsilon = epsilon
            def addPosition(self, position):
                self.position = position
            
            @property
            def getInfo(self):
                print("-"*30)
                print("name     :", self.name)
                print("element  :", self.element)
                print("index    :", self.element)
                print("residue  :", self.residue)
                print("id       :", self.id)
                print("charge   :", self.charge)
                print("sigma    :", self.sigma)
                print("epsilon  :", self.epsilon)
                print("position :", self.position)
                print("-"*30, "\n")
        
        harmonic_forces = forces[0]
        nonbonded_forces = forces[1]
        periodictorison_forces = forces[2]

        unique_force_types = []
        for force in forces:
            if type(force) not in unique_force_types:
                unique_force_types.append(type(force))
                
        print(unique_force_types)
        
        """
        [<class 'openmm.openmm.HarmonicBondForce'>, 
        <class 'openmm.openmm.NonbondedForce'>, 
        <class 'openmm.openmm.PeriodicTorsionForce'>, 
        <class 'openmm.openmm.CMMotionRemover'>,
        <class 'openmm.openmm.HarmonicAngleForce'>, 
        <class 'openmm.openmm.MonteCarloBarostat'>]
        """
        
        atoms = []
        i = 0
        for atom in topology.atoms():
            if i > 10:
                break
            else:
                atom_params = nonbonded_forces.getParticleParameters(i)
                atom_register = AtomInfo()
                atom_register.addNameAndId(atom.name, atom.element, atom.index, atom.residue, atom.id)
                atom_register.addNonbondedParams(atom_params[0], atom_params[1], atom_params[2])
                atom_register.addPosition(positions[i])
                atoms.append(atom_register)
                print(atom.name, atom.element, atom.index, atom.residue, atom.id)
                i += 1
        
        sample_0 = harmonic_forces.getBondParameters(5)
        print(sample_0)
        sample_1 = nonbonded_forces.getParticleParameters(5)
        print(sample_1)
        num_nonbonded_forced_params = nonbonded_forces.getNumParticles()
        
        for atom in atoms:
            atom.getInfo
        
        print("-"*30)
        print("Num atoms recorded in nonbonded_force :", num_nonbonded_forced_params)
        print("Num atoms recorded in topology(system):", len(atoms))
        print("Length of positions recorded in state :", len(positions))
        print("Num of force types recorded in state  :", len(forces))
        print("Length of vels recorded in state      :", len(vels))
        print("-"*30, "\n")

        """
        charge : double
            the charge of the particle, measured in units of the proton charge
        sigma : double
            the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
        epsilon : double
            the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
        """
        
        """
        Todo
        1. [ ] conenct all the information that gonna be the parts of output(interatomic-force)
            - information necessary includes
                1. [ ] atom id (the master id of each atoms registered in System)
                1. [ ] atom type(CA, H, N, etc..)
                1. [ ] atom's position (x,y,z)
                1. [ ] atom's charge
                1. [ ] atom's sigma
                1. [ ] atoms's epsilon
        1. [ ] calculate interatomic force between atoms obtained in task1 by classical for looping with cutoff.
        1. [ ] calculate interatomic force in PME method.
        1. [ ] write out the time-series date of the interatomic forces inside a molecule.
        
        """
        
        for i in range(len(atoms)):
            if i >= len(atoms):
                break
            else:
                atoms_j = atoms[i+1:]
                for j in range(len(atoms_j)):
                    self._out.write(f'{simulation_time},{atoms[i].index}_{atoms[i].name},{atoms_j[j].index}_{atoms_j[j].name},{atoms[i].charge}\n'.encode(encoding='UTF-8'))
            # self._out.write(f'{simulation_time},{},{_params[1]},{_params[2]},{_params[3]}\n'.encode(encoding='UTF-8'))

        # for bond in bonds:
        #     donor, accepter, quantity = bond[0], bond[1], 0
        #     self._out.write(f'{simulation_time},{str(donor.id)}_{donor.name},{str(accepter.id)}_{accepter.name},{quantity}\n'.encode(encoding='UTF-8'))
        
        sys.exit()
        
class InterAtomicForceReporter(InterAtomicReporter):
    
    def __init__(self) -> None:
        super().__init__()

    def describeNextReport(self, simulation: Simulation):
        return super().describeNextReport(simulation)
    
    def report(self, simulation: Simulation, state: mm.State):
        return super().report(simulation, state)
    
    