from __future__ import absolute_import
__author__ = "Yudai Okubo"
__version__ = "0.0.1"

import sys
import math

import openmm as mm
from openmm.app import Simulation, Topology
from openmm.unit import *

class AtomInfo():
    """
    Base class of the atom class for interatomic quantity calculation.
    """

    def __init__(self) -> None:
        self.name = ""
        self.element = None
        self.index = None
        self.residue = None
        self.id = None

    def add_name_and_ids(self, name, element, index, residue, id):
        self.name = name
        self.element = element
        self.index = index
        self.residue = residue
        self.id = id

    def add_nonbonded_params(self, charge, sigma, epsilon):
        self.charge = charge
        self.sigma = sigma
        self.epsilon = epsilon

    def add_position(self, position):
        self.position = position

    @property
    def get_info(self):
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
                
class Pair():
    """
    Base class of Pair 
    """

    def __init__(self, atom_i: AtomInfo, atom_j: AtomInfo, id_user_defined: int = None) -> None:
        """_summary_

        Parameters
        ----------
        atom_i : AtomInfo
            _description_
        atom_j : AtomInfo
            _description_
        id : int
            user defined id for each atoms' pairs
        """
        self.atom_i = atom_i
        self.atom_j = atom_j
        self.id = id_user_defined
    
class Pairs(object):
    """
    Pairs class
    
    # Description
    Generator class for pairwise interatomic force calculation.
    
    """
    def __init__(self) -> None:
        self.pairs: list(Pair) = []
        self._i = 0
        
    def __iter__(self):
        return self

    def __next__(self):
        if self._i== len(self.pairs):
            raise StopIteration
        pair = self.pairs[self._i]
        self._i += 1
        return pair
    
    def add_pair(self, atom_i: AtomInfo, atom_j: AtomInfo, id_user_defined: int = None) -> bool:
        try:
            self.pairs.append(Pair(atom_i, atom_j, id_user_defined))
            return True
        except:
            return False
        
    @property
    def get_num_pairs(self):
        return len(self.pairs)
    
    def setCutOffLenght(self, cutoff: nanometers):
        """setVutOffLength

        Parameters
        ----------
        cutoff : mm.unit.nanometer
            cutoff length. this value represents the threshold length between two atoms to decide whether InterAtomicForce class would calculate its interatomic force or not.
        """
        self.cutoff = cutoff
    
    def cal_interatomic_force(self):
        self._cal_interatomic_force()
    
    def _cal_interatomic_force(self) -> kilojoules/(nanometers**2*moles):
        
        for pair in self.pairs:
            
            atom_i = pair.atom_i
            atom_j = pair.atom_j
            
            charge_i = atom_i.charge
            charge_j = atom_j.charge
            sigma_i = atom_i.sigma
            sigma_j = atom_j.sigma
            epsilon_i = atom_i.sigma
            epsilon_j = atom_j.sigma
            
            pos_i = atom_i.position
            pos_j = atom_j.position
            
            r2 = (pos_i[0] - pos_j[0])**2 + (pos_i[1] - pos_j[1])**2 + (pos_i[2] - pos_j[2])**2
            r1 = r2**(1/2)
            r7 = r2**3 * r1
            r13 = r2**6 * r1
            
            sigma_ij = (sigma_i + sigma_j) / 2
            epsilon_ij = (epsilon_i*epsilon_j)**(0.5)
            
            f_lj = 4*epsilon_ij * (-12*sigma_ij**12 / r13 + 6*sigma_ij**6 / r7)*(pos_j - pos_i)/r1
            
            print(f_lj)
            f_columnb = (charge_i * charge_j / r2)*(pos_j - pos_i)/r1
            pair.fij = f_columnb + f_lj
            
    def get_interatomic_force(self, index:int ) -> kilojoules/(nanometers**2*moles):
        """get_interatomic_force

        Parameters
        ----------
        index : int
            _description_

        Returns
        -------
        inter-atomic force
            kilojoules/(nanometers**2*moles)
        """
        return self.pairs[index].fij