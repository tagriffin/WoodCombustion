# This file computes absolute enthalpies as well as average heat capacities

import numpy as np

class Molecule:
    """
    Store molecule info here
    """

    def __init__(self, name):
        """
        Pass parameters describing molecules
        """
        #! name
        self.name = name      # below are showmate coeffs from NIST, the molecular wt. and the range of validity for shomate
        self.shomate = np.array([("H2", 33.066178, -11.363417,11.432816, -2.772874, -0.158558,  -9.980797,  172.707974, 0, 2.018,"298-1000 K"),\
                   ("H2O", 30.09200, 6.832514,6.793435,  -2.534480, 0.082139, -250.8810, 223.3967, -241.8264, 18.020,"500-1700 K"), \
                   ("CO", 25.56759, 6.096130,4.054656, -2.671301, 0.131021, -118.0089, 227.3665,-110.5271, 28.010, "298-1300 K "), \
                   ("CO2low", 24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224,44.010, "298-1200 K"),\
                   ("CO2high", 58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224, 44.010, "1200-6000"),\
                   ("CH4", -0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, 158.7163, -74.87310, 16.040, "298-1300 K"),\
                   ("C2H5OH", 0, 0, 0, 0, 0, -276, 0, -276, 46.0684, "298"),\
                   ("C3H8", 0, 0, 0, 0, 0, -103.85, 0, -103.85, 44.0956, "298"),\
                   ("O2low", 31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0, 32.000, "100-700 K"),\
                   ("O2high", 30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0, 32.000, "700-2000 K"),\
                   ("N2low", 28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168, 0.0, 28.010, "100-500 K"),\
                   ("N2high", 19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0.0, 28.010, "500-2000 K")])

        row, column = np.where(self.shomate==self.name)
        self.coeffs = np.asfarray(self.shomate[row,1:10])
        self.coeffs = self.coeffs[0,:]
        self.A,self.B,self.C,self.D,self.E,self.F,self.G,self.H, self.I = self.coeffs
        self.R = 8.314 #J/(mol K)
        self.p0 = 1e5

    def h_abs(self, T): # h in kJ/mol
        """ Determining the absolute enthalpy using Shomate equation in kJ/mol"""
        t = T/1000
        return self.A*t + self.B*t**2/2 + self.C*t**3/3 + self.D*t**4/4 - self.E/t + self.F

    def h_form(self):
        return self.H

    def s_abs(self, T, p):  #pressure in Pa, s in J/mol
        t = T/1000
        return self.A*np.log(t) + self.B*t + self.C*t**2/2 + self.D*t**3/3 - self.E/(2*t**2) + self.G - self.R*np.log(p/self.p0)

    def g_abs(self, T, p):  #pressure in Pa, g in kJ/mol
        return self.h_abs(T) - T*self.s_abs(T,p)/1000

    def Cpm_ave(self,T_low,T_high): # Cpm in kJ/(kmol K)
        return (self.h_abs(T_high)-self.h_abs(T_low))/(T_high - T_low)*1000

    def Cp_ave(self,t_low,t_high):
        T_high = t_low + 273.15
        T_low = t_high + 273.15 
        return (self.h_abs(T_high)-self.h_abs(T_low))/(T_high - T_low)*1000/self.I
    
        
