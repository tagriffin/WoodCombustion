'''This file allows the calculation of various properties of the
species contained in the NASA.cti File

Timothy Griffin, FHNW Windisch, May 2021
update 21 October 2023'''

import cantera as ct
import numpy as np

ct.suppress_thermo_warnings() # prevents errors with liquid fuels from being listed

# Create a dict of species keyed by their names for entire NASA list
S = {s.name: s for s in ct.Species.list_from_file('nasa_gas.yaml')}

class NASA_Species(object):
    def __init__(self, Name):
        #! name
        self.name = Name
        self.SpecSet = []
        self.SpecSet.append(S[self.name])
        self.gas = ct.Solution(thermo='IdealGas', species = self.SpecSet)
        self.MW = self.gas.mean_molecular_weight
        self.elements = self.gas.element_names
        self.num_elements = len(self.elements)
        
        self.element_massfcn = np.zeros(self.num_elements)
        for i in range (0, self.num_elements):
            self.element_massfcn[i] = self.gas.elemental_mass_fraction(i)

        self.element_molefcn = np.zeros(self.num_elements)
        for i in range (0, self.num_elements):
            self.element_molefcn[i] = self.gas.elemental_mole_fraction(i)
            
        # Create a dictionary with elements and their mass fractions
        self.elem_mfcn = {self.elements[j]: self.element_massfcn[j] for j in range(0, self.num_elements)}

        # Create a dictionary with elements and their mole fractions
        self.elem_molfcn = {self.elements[j]: self.element_molefcn[j] for j in range(0, self.num_elements)}

    def h_mol(self, Temp, Pressure): # this returns the absolute molar enthalpy in kJ/mol
        self.gas.TP = Temp, Pressure 
        return self.gas.enthalpy_mole /1e6

    def h_mass(self, Temp, Pressure): # this returns the absolute molar enthalpy in kJ/kg
        self.gas.TP = Temp, Pressure 
        return self.gas.enthalpy_mass /1e3

    def s_mol(self, Temp, Pressure): # this returns the absolute molar entropy in J/(mol K)
        self.gas.TP = Temp, Pressure 
        return self.gas.entropy_mole /1e3

    def s_mass(self, Temp, Pressure): # this returns the absolute molar entropy in J/(mol K)
        self.gas.TP = Temp, Pressure 
        return self.gas.entropy_mass /1e3 

    def g_mol(self, Temp, Pressure): # this returns the absolute molar Gibbs energy in kJ/mol
        self.gas.TP = Temp, Pressure 
        return self.gas.gibbs_mole /1e6
    
