from __future__ import absolute_import
__author__ = "Yudai Okubo"
__version__ = "0.0.1"

import sys

import openmm as mm
from openmm.app import Simulation, Topology
import openmm.unit

class AtomInfo():
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
    
class Pairs():
    def __init__(self) -> None:
        self.pairs: list(Pair) = None
        pass
    
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
            
            atom_i = self.pair.atom_i
            atom_j = self.pair.atom_j
            
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
            
            f_lj = ()

            f_columnb = charge_i * charge_j / r2
            pair.fij = f_columnb
            
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