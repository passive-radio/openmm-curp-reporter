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