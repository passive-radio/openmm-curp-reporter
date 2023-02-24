from __future__ import absolute_import
__author__ = "Yudai Okubo"
__version__ = "0.0.1"

import sys

import openmm as mm
from openmm.app import Simulation, Topology
from openmm.unit import *

from .atominfo import AtomInfo

class InterAtomicForce():
    
    def __init__(self) -> None:
        self.cutoff = None
    
    def setCutOffLenght(self, cutoff: nanometers):
        """setVutOffLength

        Parameters
        ----------
        cutoff : mm.unit.nanometer
            cutoff length. this value represents the threshold length between two atoms to decide whether InterAtomicForce class would calculate its interatomic force or not.
        """
        self.cutoff = cutoff
        self.donor = None
        self.accepter = None
        self.force = None
        self.id = None
        
    def cal_interatomic_force(self, atom_i: AtomInfo, atom_j: AtomInfo) -> kilojoule/(nanometer**2*mole):

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
        
        return f_columnb