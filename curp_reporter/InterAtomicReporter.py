from __future__ import absolute_import
__author__ = "Yudai Okubo"
__version__ = "0.0.1"

import sys

import openmm as mm
from openmm.app import Simulation, Topology
from openmm.unit import *

from .interatomic_class import AtomInfo, Pairs, Pair

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
        """_summary_

        Parameters
        ----------
        simulation : Simulation
            _description_
        state : mm.State
            _description_

        Returns
        -------
        _type_
            _description_
        """
class InterAtomicForceReporter(InterAtomicReporter):
    
    def __init__(self, file, reportInterval) -> None:
        super().__init__(file, reportInterval)
        self._out.write('step,donor,accepter,force\n'.encode('UTF-8'))

    def describeNextReport(self, simulation: Simulation):
        return super().describeNextReport(simulation)
    
    def report(self, simulation: Simulation, state: mm.State):
        """_summary_

        Parameters
        ----------
        simulation : Simulation
            _description_
        state : mm.State
            _description_

        Returns
        -------
        _type_
            _description_
        """
        
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
        
        bonds = []
        for bond in topology.bonds():
            bonds.append(bond)
        
        unique_force_types = []
        for force in forces:
            if type(force) not in unique_force_types:
                unique_force_types.append(type(force))
                
        # print(unique_force_types)
        for force in forces:
            if type(force) == mm.NonbondedForce:
                nonbonded_forces = force
            elif type(force) == mm.HarmonicBondForce:
                harmonic_bond_forces = force
        
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
                atom_register.add_name_and_ids(atom.name, atom.element, atom.index, atom.residue, atom.id)
                atom_register.add_nonbonded_params(atom_params[0], atom_params[1], atom_params[2])
                atom_register.add_position(positions[i])
                atoms.append(atom_register)
                i += 1
                
        
        for atom in atoms:
            atom.get_info
            
        
        # print("-"*30)
        # print("Num atoms recorded in nonbonded_force :", nonbonded_forces.getNumParticles())
        # print("Num atoms recorded in topology(system):", len(atoms))
        # print("Length of positions recorded in state :", len(positions))
        # print("Num of force types recorded in state  :", len(forces))
        # print("Length of vels recorded in state      :", len(vels))
        # print("Num bonds recorded in harmonic_bond_force:", harmonic_bond_forces.getNumBonds())
        # print("Num bonds recorded in topology        :", len(bonds))
        # print("-"*30, "\n")
        
        # sample_harmonic_bond = harmonic_bond_forces.getBondParameters(1009)
        # sample_bond = bonds[1009]
        # print(sample_bond)
        # print(sample_harmonic_bond)
        

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
        1. [x] conenct all the information that gonna be the parts of output(interatomic-force)
            - information necessary includes
                1. [x] atom id (the master id of each atoms registered in System)
                1. [x] atom type(CA, H, N, etc..)
                1. [x] atom's position (x,y,z)
                1. [x] atom's charge
                1. [x] atom's sigma
                1. [x] atoms's epsilon
        1. [ ] calculate interatomic force between atoms obtained in task1 by classical for looping with cutoff.
        1. [ ] calculate interatomic force in PME method.
        1. [ ] write out the time-series date of the interatomic forces inside a molecule.
        """
        
        pairs = Pairs()
        for i in range(len(atoms)):
            if i >= len(atoms):
                break
            else:
                atoms_j = atoms[i+1:]
                for j in range(len(atoms_j)):
                    pairs.add_pair(atoms[i], atoms_j[j])
                    
        pairs.cal_interatomic_force()
        
        for pair in pairs:
            self._out.write(f'{simulation_time},{pair.atom_i.index}_{pair.atom_i.name},{pair.atom_j.index}_{pair.atom_j.name},{pair.fij}\n'.encode(encoding='UTF-8'))
        
        # for bond in bonds:
        #     donor, accepter, quantity = bond[0], bond[1], 0
        #     self._out.write(f'{simulation_time},{str(donor.id)}_{donor.name},{str(accepter.id)}_{accepter.name},{quantity}\n'.encode(encoding='UTF-8'))
        
        return super().report(simulation, state)
    
    