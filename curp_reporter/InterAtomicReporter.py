from __future__ import absolute_import
__author__ = "Yudai Okubo"
__version__ = "0.0.1"

import openmm as mm
from openmm.app import Simulation, Topology
import openmm.unit

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
            
        self._out.write('step,donor,accepter,quantity\n'.encode('UTF-8'))
        
    def describeNextReport(self, simulation: Simulation):
        steps = self.__reportInterval - simulation.currentStep%self.__reportInterval
        return (steps, self, True, False, False, False)
    
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
        positions = state.getPositions()
        
        system = simulation.system
        simulation_time = state.getTime()
        topology = simulation.topology
        bonds = topology.bonds()
        
        numForces = system.getNumForces()
        forces = system.getForces()
        
        for i in topology.atoms():
            print(i)
        
        harmonic_forces = forces[0]
        num_harmonic_forces = harmonic_forces.getNumBonds()
        
        harmonic_forces_params = []
        for i in range(num_harmonic_forces):
            _params = harmonic_forces.getBondParameters(i)
            self._out.write(f'{simulation_time},{_params[0]},{_params[1]},{_params[2]},{_params[3]}\n'.encode(encoding='UTF-8'))


        # for bond in bonds:
        #     donor, accepter, quantity = bond[0], bond[1], 0
        #     self._out.write(f'{simulation_time},{str(donor.id)}_{donor.name},{str(accepter.id)}_{accepter.name},{quantity}\n'.encode(encoding='UTF-8'))
        
        
        
class InterAtomicForceReporter(InterAtomicReporter):
    
    def __init__(self) -> None:
        super().__init__()

    def describeNextReport(self, simulation: Simulation):
        return super().describeNextReport(simulation)
    
    def report(self, simulation: Simulation, state: mm.State):
        return super().report(simulation, state)
    
    