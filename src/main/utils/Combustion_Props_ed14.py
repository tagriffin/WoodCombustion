'''This file creates an object Fuel_1 or Fuel_2 and has a number of instances relevant for
combustion of Type 1 and Type 2 Fuels at 1 bar total pressure.

For Type 1 Fuels, a mixture of fuels can be used, with mole fractions adding up to 1
The names must match those used in the Cantera NASA.cti file. 

For Type 2 Fuels water to be added,
for example to reduce NOx in oil combustion, or to accommodate the water content in biomass fuels

Timothy Griffin, FHNW Windisch, May 2021 ... Updated 29 November 2021'''

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd
from pandas import ExcelWriter

import cantera as ct
from src.main.utils import Properties_NASA_v5 as NASA
from src.main.utils import A8
import CoolProp.CoolProp as CP

H2 = NASA.NASA_Species('H2')
H2O = NASA.NASA_Species('H2O')
SO2 = NASA.NASA_Species('SO2')
N2 = NASA.NASA_Species('N2')
O2 = NASA.NASA_Species('O2')
CO2 = NASA.NASA_Species('CO2')

zeta_O2_Air = 0.2314
zeta_N2_Air = 1 - zeta_O2_Air

type = {'CH4': 'Typ 1', 'C2H6': 'Typ 1', 'C3H8': 'Typ 1', 'C4H10': 'Typ 1', 'H2': 'Typ 1',
         'Benzin': 'Typ 2', 'Diesel': 'Typ 2', 'Holz': 'Typ 2', 'Erdgas': 'Typ 1'}

def Fuel_Type(Name):
    return type[Name]

# Class for pure or mixture fuels, type 1
class Fuel_1(object):
    def __init__(self, name, Components, MolFcns): # NASA_name, Vector(Components, Vector (molfcs)
        """
        Pass parameters describing molecules
        """
        self.name = name
        # Cantera must have a "mixture" here: if just one component must "add" N2
        if len(Components) == 1:
            Components.append('N2')
            MolFcns.append(0.0)
        self.comp = Components
        self.molfcn = MolFcns

        # Setting the fuel mixture in Cantera
        # Create a dict of species keyed by their names for entire NASA list
        S = {s.name: s for s in ct.Species.listFromFile('nasa.xml')}
        SpecSet = []
        # Adding the components for the object fuel
        for i in range (len(self.comp)):
            SpecSet.append(S[self.comp[i]])
                           
        Fuel = ct.Solution(thermo='IdealGas', species = SpecSet)
        Fuel.X = self.molfcn
        self.MixFuel = Fuel
        
        # extracting the masses of the elements in the fuel
        self.C = 0; self.Cm = 0
        self.H = 0; self.Hm = 0
        self.O = 0; self.Om = 0
        self.N = 0; self.Nm = 0
        self.S = 0; self.Sm = 0

        self.mw = Fuel.mean_molecular_weight # Molecular wt in kg/kmol
        
        for i in range (len(Fuel.element_names)):
            if Fuel.element_names[i] == 'C':
                self.C = Fuel.elemental_mass_fraction(i)
                self.Cm = Fuel.elemental_mole_fraction(i)
            if Fuel.element_names[i] == 'H':
                self.H = Fuel.elemental_mass_fraction(i)
                self.Hm = Fuel.elemental_mole_fraction(i)
            if Fuel.element_names[i] == 'O':
                self.O = Fuel.elemental_mass_fraction(i)
                self.Om = Fuel.elemental_mole_fraction(i)
            if Fuel.element_names[i] == 'N':
                self.N = Fuel.elemental_mass_fraction(i)
                self.Nm = Fuel.elemental_mole_fraction(i)
            if Fuel.element_names[i] == 'S':
                self.S = Fuel.elemental_mass_fraction(i)
                self.Sm = Fuel.elemental_mole_fraction(i)

        # MW of Fuel with C1
        self.mw_C1 = self.Cm*12.0107 + self.Hm*1.00794 + self.Om*15.994 + \
                  self.Nm*14.0067 + self.Sm*32.0
        self.C_Factor = self.mw/self.mw_C1

        # Coefficients for balancing the molar reaction for complete oxidation of wet fuel
        # C(Cm)H(Hm)S(Sm)O(Om)N(Nm) + lambda*a(O2 + 3.76N2) = bCO2 + cH2O + dSO2 + (lambda*e + f)N2 + (lambda-1)aO2
        self.a = (2*self.Cm + self.Hm/2 + 2*self.Sm - self.Om)/2
        self.b = self.Cm 
        self.c = self.Hm/2 
        self.d = self.Sm 
        self.e = self.a * 3.76 
        self.f = self.Nm/2 
        
        print('The fuel is C_%0.3f' %(self.Cm),'H_%0.3f' %(self.Hm), 'O_%0.3f' %(self.Om))

    def Hu(self): # in MJ/kg  fuel
        # for fuels without S and N
        T0 = 298.15 # ref Temp in K
        p0 = 1e5 # ref Pressure in Pa
        Fuel = self.MixFuel
        Fuel.TP = T0, p0
        Hu_molar = Fuel.enthalpy_mole /1e6 - self.C_Factor * (self.b*CO2.h_mol(T0,p0) + self.c*H2O.h_mol(T0,p0) + self.d*SO2.h_mol(T0,p0))
        return Hu_molar/self.mw

    def Ho(self): # higher heating value in MJ/kg fuel at 25 °C(Brennwert)
        mu_H2O = self.c * self.C_Factor * 18.02/self.mw
        return self.Hu() + mu_H2O * (2.4423)

    def Omin(self): # this is calculated for 1 kmol fuel
        return self.C_Factor * (self.Cm + 1/4 * self.Hm + self.Sm - 0.5*self.Om)

    def omin(self): # this is calculated for 1 kg fuel
        return self.Omin() * 31.9988/self.mw

    def lmin(self): # this is calculated for 1 kg fuel
        return self.omin()/zeta_O2_Air

    def lmin_wet(self, Tair_in, Rel_humidity): # this is calculated for 1 kg fuel and moist air and 1 bar pressure
        p = 1 # 1 bar combustion pressure
        Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
        x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        return (self.omin()/zeta_O2_Air) * (1 + x_humidity)

    def Lmin(self): # this is calculated for 1 kmol wet fuel
        return self.Omin() / 0.2095

    # Exhaust Gas composition (only for Lambda > 1), using the molar comp of fuel
    def x_CO2(self, Lambda):
        return self.b/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_CO2(self, Lambda): # works only for fuel without S and N!
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total = m_CO2 + m_H2O + m_O2ex + m_N2
        return m_CO2/mass_total

    def x_CO2_tr(self, Lambda):
        return self.b/(self.b + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def x_H2O(self, Lambda):
        return self.c/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_H2O(self, Lambda):
        m_Air = self.lmin() * Lambda 
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total =  m_CO2 + m_H2O + m_O2ex + m_N2
        return m_H2O/mass_total

    def TauPt(self, Lambda, pressure): # pressure in Pa, only works for combustion in dry air
        p_w = self.x_H2O(Lambda) * pressure
        return CP.PropsSI("T","P", p_w,"Q",1,"water")

    def x_SO2(self, Lambda):
        return self.d/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def x_N2(self, Lambda):
        return (Lambda*self.e + self.f)/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_N2(self, Lambda):
        m_Air = self.lmin() * Lambda 
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total =  m_CO2 + m_H2O + m_O2ex + m_N2
        return m_N2/mass_total

    def x_O2(self, Lambda):
        return ((Lambda - 1)*self.a)/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_O2(self, Lambda):
        m_Air = self.lmin() * Lambda 
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total =  m_CO2 + m_H2O + m_O2ex + m_N2
        return m_O2/mass_total

    def x_O2_tr(self, Lambda):
        return ((Lambda - 1)*self.a)/(self.b + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)
    
    def Lamb_CO2_tr(self, CO2_tr):
        def Lamb_fcn(Lamb):
            return self.x_CO2_tr(Lamb) - CO2_tr
        Lamb = 1.5 # intitial guess
        return optimize.fsolve(Lamb_fcn, Lamb)

    def Lamb_O2_tr(self, O2_tr):
        def Lamb_fcn(Lamb):
            return self.x_O2_tr(Lamb) - O2_tr
        Lamb = 1.5 # intitial guess
        return optimize.fsolve(Lamb_fcn, Lamb)

    def CO2_energy(self): # kg CO2 per MJ energy
        return self.b * CO2.MW/self.mw_C1/self.Hu()

    def CO2_kg_fuel(self): # kg CO2 per kg fuel
        return self.b * CO2.MW/self.mw_C1

    def H_in(self, Tair_in, Lambda, Rel_humidity = 0):
        # Using the Cpave data from Appendix 8
        # Problems with humidity if Tair_in > 647 K
        p = 1 # 1 bar combustion pressure; for ideal gas mixture pressure just important for humidity
        p_Pa = p * 1e5
        Fuel = self.MixFuel
        Fuel.TP = 298.15, p_Pa # assuming fuel is not preheated and equal to 25°C
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0

        #Tair_in = float(input('Enter the inlet temperature of  air in °C  ')) + 273.15
        
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet fuel
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_H2O_Air = m_Air * x_humidity
        return Fuel.enthalpy_mass/1000 + m_O2*O2.h_mass(Tair_in, p_Pa) +  \
               m_N2*N2.h_mass(Tair_in, p_Pa) + m_H2O_Air * H2O.h_mass(Tair_in, p_Pa)

    def H_out(self,Tair_in, Tout, Lambda, Rel_humidity = 0): # H_out works only  for Lambda > 1
        # Basis is 1 kg of wet fuel
        p = 1 # 1 bar combustion pressure
        p_Pa = p * 1e5
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_H2O_Air = m_Air * x_humidity 

        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H + m_H2O_Air # 
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        return m_CO2*CO2.h_mass(Tout, p_Pa) + m_H2O*H2O.h_mass(Tout, p_Pa) +  \
               m_O2ex*O2.h_mass(Tout, p_Pa) + m_N2*N2.h_mass(Tout, p_Pa)

    def H_in_A8 (self, Tair_in, Lambda, Rel_humidity = 0):
        # Assumption: can use N2* for all N2 (including the N in the fuel)
        # Problems with humidity if Tair_in > 647 K
        p = 1 # 1 bar combustion pressure; for ideal gas mixture pressure just important for humidity
        p_Pa = p * 1e5
        t_in = Tair_in - 273.15
        Fuel = self.MixFuel
        Fuel.TP = Tair_in, p_Pa
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0

        #Tair_in = float(input('Enter the inlet temperature of  air in °C  ')) + 273.15
        
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet fuel
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_H2O_Air = m_Air * x_humidity
        return self.Hu()* 1000 + (m_O2*A8.Cp_ave_O2(t_in) +  \
               m_N2*A8.Cp_ave_N2s(t_in) + m_H2O_Air * A8.Cp_ave_H2O(t_in))* (t_in - 25)

    def H_out_A8(self,Tair_in, Tout, Lambda, Rel_humidity = 0): # H_out works only  for Lambda > 1
        # Using the Cpave data from Appendix 8: only valid up to 2250 °C !
        # Basis is 1 kg of wet fuel
        # Assumption: can use N2* for all N2 (including the N in the fuel)
        p = 1 # 1 bar combustion pressure
        p_Pa = p * 1e5
        t_out = Tout - 273.15
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_H2O_Air = m_Air * x_humidity 

        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H + m_H2O_Air # 
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        return (m_CO2*A8.Cp_ave_CO2(t_out) + m_H2O*A8.Cp_ave_H2O(t_out) +  \
               m_O2ex*A8.Cp_ave_O2(t_out)+ m_N2*A8.Cp_ave_N2(t_out)) * (t_out - 25)

    def T_ad(self, Tair_in, Lambda, Rel_humidity = 0):  # Adiabatic flame temperature calculation (full oxidation)
        # Basis is 1 kg of wet fuel
        def H_fcn(T):
            return self.H_out(Tair_in, T, Lambda, Rel_humidity)- self.H_in(Tair_in, Lambda, Rel_humidity)
        T = 1200 # Initial Guess for Tad in K
        return optimize.fsolve(H_fcn, T)

    def T_ad_A8(self, Tair_in, Lambda, Rel_humidity = 0):  # Adiabatic flame temperature calculation (full oxidation)
        # Basis is 1 kg of wet fuel
        # Using the Cpave data from Appendix 8: only valid up to 2250 °C !
        def H_fcn(T):
            return self.H_out_A8(Tair_in, T, Lambda, Rel_humidity)- self.H_in_A8(Tair_in, Lambda, Rel_humidity)
        T = 1200 # Initial Guess for Tad in K
        return optimize.fsolve(H_fcn, T)
    
    def Cp_ave_stExhaust(self,t_low,t_high):  # this is calculated for dry fuel (water content = 0)
        p_Pa = 1e5
        T_high = t_low + 273.15
        T_low = t_high + 273.15
        
        # masses of stoichiometric exhaust gas
        m_Air = self.lmin() # mass of stoichiometric dry air per kg wet wood 
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H 
        m_N2 = m_Air * zeta_N2_Air
        m_total = m_Air + 1 # mass of stoichiometric exhaust gas of 1 kg of wet wood
            
        h_high = m_CO2*CO2.h_mass(T_high, p_Pa) + m_H2O*H2O.h_mass(T_high, p_Pa) + m_N2*N2.h_mass(T_high, p_Pa)
        h_low = m_CO2*CO2.h_mass(T_low, p_Pa) + m_H2O*H2O.h_mass(T_low, p_Pa) + m_N2*N2.h_mass(T_low, p_Pa)
        return (h_high -h_low)/(T_high - T_low)/m_total

    def Cp_ave_plot(self):
        temps = np.linspace(0,2200,200) # the range of data is 0 °C  to 2'200 °C (the A8 data go to -60 °C)
        Cp_ave = np.zeros(len(temps))
        for i in range(0,len(temps)):
            Cp_ave[i] = self.Cp_ave_stExhaust(25,temps[i])
        plt.plot(temps,Cp_ave, label = 'Cp_ave_stoich_exhaust')
        plt.title('Cp_average from 25 °C to t(°C) for stoich exhaust for ' + self.name)
        plt.xlabel('Final Temperature [deg C]')
        plt.ylabel('Cp_ave in kJ/(kg K)')
        plt.legend(loc='best')
        plt.grid()
        plt.show()

    def Write_Cp_ave(self):
        temps = np.linspace(0,2200,221) # the range of data is 0 °C  to 2'200 °C
        Cp_ave = np.zeros(len(temps))
        for i in range(0,len(temps)):
            Cp_ave[i] = self.Cp_ave_stExhaust(25,temps[i])
        writer = pd.ExcelWriter(self.name + '_Cp_ave_stExhaust.xlsx', engine='xlsxwriter')
        df1 = pd.DataFrame({'Temp (°C)': temps})
        df1.to_excel(writer, sheet_name=self.name, index=False)
        df2 = pd.DataFrame({'Cp_ave_stExhaust': Cp_ave})
        df2.to_excel(writer, sheet_name=self.name, startcol=1, index=False)
        writer.save()

class Fuel_2(object):
    def __init__(self, name, gamma, Hu_dry = 0, water_content = 0): # mass fractions of elements, Hu per kg dry fuel
        """
        Pass parameters describing molecules
        """
        #! name
        self.name = name
        # extracting the masses of the elements in the dry fuel
        self.C = 0
        self.H = 0
        self.O = 0
        self.N = 0
        self.S = 0
        for i in range (len(gamma)):
            if gamma[i] == 'C':
                self.C = gamma[i+1]
            if gamma[i] == 'H':
                self.H = gamma[i+1]
            if gamma[i] == 'O':
                self.O = gamma[i+1]
            if gamma[i] == 'N':
                self.N = gamma[i+1]
            if gamma[i] == 'S':
                self.S = gamma[i+1]
                
        # normalizing these values to obtain mass fractions of dry fuel
        ElTotal = self.C + self.H + self.O + self.N + self.S
        self.C = self.C/ElTotal
        self.H = self.H/ElTotal
        self.O = self.O/ElTotal
        self.N = self.N/ElTotal
        self.S = self.S/ElTotal

        self.w = water_content
        self.u = self.w/(1 - self.w) # 1 kg of dry wood has u kg water
        self.Hu_dry = Hu_dry

        # normalized mass fractions of wet fuel
        Total_w = self.C + self.H + self.O + self.N + self.S + self.u
        self.C_w = self.C/Total_w
        self.H_w = (self.H + self.u *(2*1.00794/18.0149)) /Total_w
        self.O_w = (self.O + self.u *(15.994/18.0149)) /Total_w
        self.N_w = self.N/Total_w
        self.S_w = self.S/Total_w
        self.Total_w = Total_w

        # molar composition of wet wood, basis 1 kg dry fuel 
        Cm_w = self.C/12.0107
        Hm_w = (self.H + self.u * (2*1.00794/18.0149))/1.00794
        Om_w = (self.O + self.u * (15.994/18.0149))/15.994
        Nm_w = self.N/14.0067
        Sm_w = self.S/32.0
        ElTotalm = Cm_w + Hm_w + Om_w + Nm_w + Sm_w
        
        # Molar composition of wet fuel, basis is C1; 
        self.Cm_w_norm = Cm_w/ElTotalm
        self.Cm_w = 1.0
        self.Hm_w = (Hm_w/ElTotalm)/self.Cm_w_norm
        self.Om_w = (Om_w/ElTotalm)/self.Cm_w_norm
        self.Nm_w = (Nm_w/ElTotalm)/self.Cm_w_norm
        self.Sm_w = (Sm_w/ElTotalm)/self.Cm_w_norm
        self.mw = self.Cm_w*12.0107 + self.Hm_w*1.00794 + self.Om_w*15.994 + \
                  self.Nm_w*14.0067 + self.Sm_w*32.0

        # Coefficients for balancing the molar reaction for complete oxidation of wet fuel
        # C(Cm)H(Hm)S(Sm)O(Om)N(Nm) + lambda*a(O2 + 3.76N2) = bCO2 + cH2O + dSO2 + (lambda*e + f)N2 + (lambda-1)aO2
        self.a = (2*self.Cm_w + self.Hm_w/2 + 2*self.Sm_w - self.Om_w)/2
        self.b = self.Cm_w
        self.c = self.Hm_w/2
        self.d = self.Sm_w
        self.e = self.a * 3.76
        self.f = self.Nm_w/2
        
        print('The fuel is C_%0.3f' %(self.Cm_w),'H_%0.3f' %(self.Hm_w), 'O_%0.3f' %(self.Om_w), \
              'with "Wassergehalt" of %0.2f' %(self.w))
        
    def omin(self): # this is calculated for 1 kg wet fuel
        return 2.664*self.C_w + 7.937*(self.H/self.Total_w) + 0.998*self.S_w - (self.O/self.Total_w)

    def Omin(self): # this is calculated for 1 kmol wet fuel
        return self.omin()* self.mw/31.988

    def lmin(self): # this is calculated for 1 kg wet fuel
        return self.omin()/zeta_O2_Air

    def lmin_dry(self): # this calculates lmin for 1 kg dry wood (as in Thermo script)
        return self.lmin()* (1 + self.u)

    def lmin_wet(self, Tair_in, Rel_humidity): # this is calculated for 1 kg fuel and moist air and 1 bar pressure
        p = 1 # 1 bar combustion pressure
        Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
        x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        return (self.omin()/zeta_O2_Air) * (1 + x_humidity)

    def Lmin(self): # this is calculated for 1 kmol wet fuel
        return self.lmin() * self.mw/29.871

    def Hu_nass(self): # in MJ/kg wet fuel
        return self.Hu_dry*(1-self.w) - 2.443*self.w

    def Hf_nass(self): # in MJ/kg wet fuel
        # for fuels without S and N
        T0 = 298.15 # ref Temp in K
        p0 = 1e5 # ref Pressure in Pa
        Hf_nassm = self.Cm_w*CO2.h_mol(T0,p0) + (self.Hm_w)/2*H2O.h_mol(T0,p0) + self.Hu_nass()*self.mw 
        return Hf_nassm / self.mw
    
    def Ho_nass(self): # higher heating value in MJ/kg fuel at 25 °C(Brennwert)
        mu_H2O = self.c * 18.02/self.mw
        return self.Hu_nass() + mu_H2O * (2.4423)
    
        # Exhaust Gas composition (only for Lambda > 1), using the molar comp of wet wood
    def x_CO2(self, Lambda):
        return self.b/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_CO2(self, Lambda): # works only for fuel without S and N!
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total =  m_CO2 + m_H2O + m_O2ex + m_N2
        return m_CO2/mass_total

    def x_CO2_tr(self, Lambda):
        return self.b/(self.b + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def x_H2O(self, Lambda):
        return self.c/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_H2O(self, Lambda): # works only for fuel without S and N!
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total =  m_CO2 + m_H2O + m_O2ex + m_N2
        return m_H2O/mass_total

    def TauPt(self, Lambda, pressure): # pressure in Pa, only works for combustion in dry air
        p_w = self.x_H2O(Lambda) * pressure
        return CP.PropsSI("T","P", p_w,"Q",1,"water")

    def x_SO2(self, Lambda):
        return self.d/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def x_N2(self, Lambda):
        return (Lambda*self.e + self.f)/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_N2(self, Lambda): # works only for fuel without S and N!
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total =  m_CO2 + m_H2O + m_O2ex + m_N2
        return m_N2/mass_total

    def x_O2(self, Lambda):
        return ((Lambda - 1)*self.a)/(self.b + self.c + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def massfcn_O2(self, Lambda): # works only for fuel without S and N!
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H
        m_O2 = m_Air * zeta_O2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        mass_total =  m_CO2 + m_H2O + m_O2ex + m_N2
        return m_O2ex/mass_total


    def x_O2_tr(self, Lambda):
        return ((Lambda - 1)*self.a)/(self.b + self.d + (Lambda*self.e + self.f) + (Lambda-1)*self.a)

    def Lamb_CO2_tr(self, CO2_tr):
        def Lamb_fcn(Lamb):
            return self.x_CO2_tr(Lamb) - CO2_tr
        Lamb = 1.5 # intitial guess
        return optimize.fsolve(Lamb_fcn, Lamb)

    def Lamb_O2_tr(self, O2_tr):
        def Lamb_fcn(Lamb):
            return self.x_O2_tr(Lamb) - O2_tr
        Lamb = 1.5 # intitial guess
        return optimize.fsolve(Lamb_fcn, Lamb)

    def CO2_energy(self): # kg CO2 per MJ energy
        return self.b * CO2.MW/self.mw/self.Hu_nass()

    def CO2_kg_fuel(self): # kg CO2 per kg fuel
        return self.b * CO2.MW/self.mw
    
    def H_in (self, Tair_in, Lambda, Rel_humidity = 0):
        p = 1 # 1 bar combustion pressure
        p_Pa = p * 1e5
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0

        #Tair_in = float(input('Enter the inlet temperature of  air in °C  ')) + 273.15
        
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet fuel
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_H2O_Air = m_Air * x_humidity
        # assumes that fuel is not preheated, thus can use Hf at same temperature as heating value
        return self.Hf_nass()*1000 + m_O2*O2.h_mass(Tair_in, p_Pa) +  \
               m_N2*N2.h_mass(Tair_in, p_Pa) + m_H2O_Air * H2O.h_mass(Tair_in, p_Pa)

    def H_out(self,Tair_in, Tout, Lambda, Rel_humidity = 0): # H_out works only  for Lambda > 1
        # Basis is 1 kg of wet fuel
        p = 1 # 1 bar combustion pressure
        p_Pa = p * 1e5
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_H2O_Air = m_Air * x_humidity 

        m_CO2 = 3.664 * self.C_w
        m_H2O = 8.937 * (self.H/self.Total_w) + self.w + m_H2O_Air # must split the H // mass water from 1kg wet wood
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        return m_CO2*CO2.h_mass(Tout, p_Pa) + m_H2O*H2O.h_mass(Tout, p_Pa) +  \
               m_O2ex*O2.h_mass(Tout, p_Pa) + m_N2*N2.h_mass(Tout, p_Pa)

    def H_in_A8(self, Tair_in, Lambda, Rel_humidity = 0):
        # Assumption: can use N2* for all N2 (including the N in the fuel)
        p = 1 # 1 bar combustion pressure
        p_Pa = p * 1e5
        t_in = Tair_in - 273.15
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0

        #Tair_in = float(input('Enter the inlet temperature of  air in °C  ')) + 273.15
        
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet fuel
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_H2O_Air = m_Air * x_humidity
        # assumes that fuel is not preheated, thus can use Hf at same temperature as heating value
        return self.Hu_nass()* 1000 + (m_O2*A8.Cp_ave_O2(t_in) +  \
               m_N2*A8.Cp_ave_N2s(t_in) + m_H2O_Air * A8.Cp_ave_H2O(t_in))* (t_in - 25)

    def H_out_A8(self,Tair_in, Tout, Lambda, Rel_humidity = 0): # H_out works only  for Lambda > 1
        # Basis is 1 kg of wet fuel
        # Assumption: can use N2* for all N2 (including the N in the fuel)
        p = 1 # 1 bar combustion pressure
        p_Pa = p * 1e5
        t_out = Tout - 273.15
        if Tair_in < 647:
            Psat = CP.PropsSI("P","T",Tair_in,"Q",1,"water")/10**5  # in bar
            x_humidity = 0.622*Rel_humidity*Psat/(p-Rel_humidity*Psat)
        else:
            x_humidity = 0
        m_Air = self.lmin() * Lambda # mass of dry air per kg wet wood
        m_H2O_Air = m_Air * x_humidity 

        m_CO2 = 3.664 * self.C_w
        m_H2O = 8.937 * (self.H/self.Total_w) + self.w + m_H2O_Air # must split the H // mass water from 1kg wet wood
        m_O2 = m_Air * zeta_O2_Air
        m_N2 = m_Air * zeta_N2_Air
        m_O2ex = m_O2 - self.lmin()* zeta_O2_Air
        return (m_CO2*A8.Cp_ave_CO2(t_out) + m_H2O*A8.Cp_ave_H2O(t_out) +  \
               m_O2ex*A8.Cp_ave_O2(t_out) + m_N2*A8.Cp_ave_N2(t_out)) * (t_out - 25)

    def T_ad(self, Tair_in, Lambda, Rel_humidity = 0):  # Adiabatic flame temperature calculation (full oxidation)
        # Basis is 1 kg of wet fuel
        def H_fcn(T):
            return self.H_out(Tair_in, T, Lambda, Rel_humidity)- self.H_in(Tair_in, Lambda, Rel_humidity)
        T = 1200 # Initial Guess for Tad in K
        return optimize.fsolve(H_fcn, T)

    def T_ad_A8(self, Tair_in, Lambda, Rel_humidity = 0):  # Adiabatic flame temperature calculation (full oxidation)
        # Basis is 1 kg of wet fuel
        # Using the Cpave data from Appendix 8: only valid up to 2250 °C !
        def H_fcn(T):
            return self.H_out_A8(Tair_in, T, Lambda, Rel_humidity)- self.H_in_A8(Tair_in, Lambda, Rel_humidity)
        T = 1200 # Initial Guess for Tad in K
        return optimize.fsolve(H_fcn, T)

    def Cp_ave_stExhaust(self,t_low,t_high):  # this is calculated for dry fuel (water content = 0)
        p_Pa = 1e5
        T_high = t_low + 273.15
        T_low = t_high + 273.15
        
        # masses of stoichiometric exhaust gas
        m_Air = self.lmin() # mass of stoichiometric dry air per kg wet wood 
        m_CO2 = 3.664 * self.C
        m_H2O = 8.937 * self.H 
        m_N2 = m_Air * zeta_N2_Air
        m_total = m_Air + 1 # mass of stoichiometric exhaust gas of 1 kg of wet wood
            
        h_high = m_CO2*CO2.h_mass(T_high, p_Pa) + m_H2O*H2O.h_mass(T_high, p_Pa) + m_N2*N2.h_mass(T_high, p_Pa)
        h_low = m_CO2*CO2.h_mass(T_low, p_Pa) + m_H2O*H2O.h_mass(T_low, p_Pa) + m_N2*N2.h_mass(T_low, p_Pa)
        return (h_high -h_low)/(T_high - T_low)/m_total

    def Cp_ave_plot(self):
        temps = np.linspace(0,2200,221) # the range of data is 0 °C  to 2'200 °C (the A8 data go to -60 °C)
        Cp_ave = np.zeros(len(temps))
        for i in range(0,len(temps)):
            Cp_ave[i] = self.Cp_ave_stExhaust(25,temps[i])
        plt.plot(temps,Cp_ave, label = 'Cp_ave_stoich_exhaust')
        plt.title('Cp_average from 25 °C to t(°C) for stoich exhaust for ' + self.name)
        plt.xlabel('Final Temperature [deg C]')
        plt.ylabel('Cp_ave in kJ/(kg K)')
        plt.legend(loc='best')
        plt.grid()
        plt.show()

    def Write_Cp_ave(self):
        temps = np.linspace(0,2200,221) # the range of data is 0 °C  to 2'200 °C
        Cp_ave = np.zeros(len(temps))
        for i in range(0,len(temps)):
            Cp_ave[i] = self.Cp_ave_stExhaust(25,temps[i])
        writer = pd.ExcelWriter(self.name + '_Cp_ave_stExhaust.xlsx', engine='xlsxwriter')
        df1 = pd.DataFrame({'Temp (°C)': temps})
        df1.to_excel(writer, sheet_name=self.name, index=False)
        df2 = pd.DataFrame({'Cp_ave_stExhaust': Cp_ave})
        df2.to_excel(writer, sheet_name=self.name, startcol=1, index=False)
        writer.save()


    
