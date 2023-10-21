''' This calculates the average mass-specific heat capacities of main
combustion exhaust gases'''

from src.main.utils import IG_Props_Comb as Props
import numpy as np
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

df1 = pd.read_excel('Thd Skript Anhang Tabelle A-08 Rev07.xlsx', sheet_name='A-8.1')
t1 = np.array(df1['T (°C)'])
Cp_Air = np.array(df1['Air'])
Cp_N2s = np.array(df1['N2*'])
Cp_N2 = np.array(df1['N2'])
Cp_O2 = np.array(df1['O2'])
Cp_CO2 = np.array(df1['CO2'])
Cp_H2O = np.array(df1['H2O'])
Cp_SO2 = np.array(df1['SO2'])

# Interpolating data from the data in table
def Cp_ave_Air(t_value):
    Cp_value = interp1d(t1,Cp_Air)
    return Cp_value(t_value)

def Cp_ave_N2s(t_value):
    Cp_value = interp1d(t1,Cp_N2s)
    return Cp_value(t_value)

def Cp_ave_N2(t_value):
    Cp_value = interp1d(t1,Cp_N2)
    return Cp_value(t_value)

def Cp_ave_O2(t_value):
    Cp_value = interp1d(t1,Cp_O2)
    return Cp_value(t_value)

def Cp_ave_CO2(t_value):
    Cp_value = interp1d(t1,Cp_CO2)
    return Cp_value(t_value)

def Cp_ave_H2O(t_value):
    Cp_value = interp1d(t1,Cp_H2O)
    return Cp_value(t_value)

def Cp_ave_SO2(t_value):
    Cp_value = interp1d(t1,Cp_SO2)
    return Cp_value(t_value)

df2 = pd.read_excel('Thd Skript Anhang Tabelle A-08 Rev07.xlsx', sheet_name='A-8.2')
t2 = np.array(df2['T (°C)'])
s_Air = np.array(df2['Air'])
s_N2s = np.array(df2['N2*'])
s_N2 = np.array(df2['N2'])
s_O2 = np.array(df2['O2'])
s_CO2 = np.array(df2['CO2'])
s_H2O = np.array(df2['H2O'])
s_SO2 = np.array(df2['SO2'])

# Interpolating data from the data in table
def s_abs_Air(t_value):
    s_abs_value = interp1d(t2,s_Air)
    return s_abs_value(t_value)

def s_abs_N2s(t_value):
    s_abs_value = interp1d(t2,s_N2s)
    return s_abs_value(t_value)

def s_abs_N2(t_value):
    s_abs_value = interp1d(t2,s_N2)
    return s_abs_value(t_value)

def s_abs_O2(t_value):
    s_abs_value = interp1d(t2,s_O2)
    return s_abs_value(t_value)

def s_abs_CO2(t_value):
    s_abs_value = interp1d(t2,s_CO2)
    return s_abs_value(t_value)

def s_abs_H2O(t_value):
    s_abs_value = interp1d(t2,s_H2O)
    return s_abs_value(t_value)

def s_abs_SO2(t_value):
    s_abs_value = interp1d(t2,s_SO2)
    return s_abs_value(t_value)

# Calculating values using Shomate polynomials(c-Values)
# Data for Gases from IG Props

def Cp_ave_Air_c(t_low,t_high):
    return 0.7686 * Cp_ave_N2_c(t_low,t_high) + 0.2314 * Cp_ave_O2_c(t_low,t_high)

def Cp_ave_N2_c(t_low,t_high):
    T_high = t_low + 273.15
    T_low = t_high + 273.15
    if T_low < 500:
        N2 = Props.Molecule('N2low') #valid 100 to 500 K
        h_low = N2.h_abs(T_low)
    else: 
        N2 = Props.Molecule('N2high') #valid 500 to 2000 K
        h_low = N2.h_abs(T_low)
    if T_high < 500:
        N2 = Props.Molecule('N2low') #valid 100 to 500 K
        h_high = N2.h_abs(T_high)
    else: 
        N2 = Props.Molecule('N2high') #valid 500 to 2000 K
        h_high = N2.h_abs(T_high)        
    return (h_high- h_low)/(T_high - T_low)*1000/N2.I

def Cp_ave_O2_c(t_low,t_high):
    T_high = t_low + 273.15
    T_low = t_high + 273.15
    if T_low < 700:
        O2 = Props.Molecule('O2low') #valid 100 to 700 K
        h_low = O2.h_abs(T_low)
    else: 
        O2 = Props.Molecule('O2high') #valid 1200 to 6000 K
        h_low = O2.h_abs(T_low)
    if T_high < 700:
        O2 = Props.Molecule('O2low') #valid 100 to 700 K
        h_high = O2.h_abs(T_high)
    else: 
        O2 = Props.Molecule('O2high') #valid 1200 to 6000 K
        h_high = O2.h_abs(T_high)
    return (h_high- h_low)/(T_high - T_low)*1000/O2.I

def Cp_ave_CO2_c(t_low,t_high):
    T_high = t_low + 273.15
    T_low = t_high + 273.15
    if T_low  < 1200:
        CO2 = Props.Molecule('CO2low') #valid 298 to 1200 K
        h_low = CO2.h_abs(T_low)
    else: 
        CO2 = Props.Molecule('CO2high') #valid 1200 to 6000 K
        h_low = CO2.h_abs(T_low)
    if T_high  < 1200:
        CO2 = Props.Molecule('CO2low') #valid 298 to 1200 K
        h_high = CO2.h_abs(T_high)
    else: 
        CO2 = Props.Molecule('CO2high') #valid 1200 to 6000 K
        h_high = CO2.h_abs(T_high)                     
    return (h_high- h_low)/(T_high - T_low)*1000/CO2.I

def Cp_ave_H2O_c(t_low,t_high):
    T_high = t_low + 273.15
    T_low = t_high + 273.15
    H2O = Props.Molecule('H2O') # valid 500 to 1700 K
    return H2O.Cp_ave(t_low,t_high)

# Plotting the data and comparing (starting always from 0 °C)
temps = np.linspace(-60,2200,200) # the range of data is -60 °C  to 2'200 °C
Cp_ave_Air_cal = np.zeros(len(temps))
Cp_ave_N2_cal = np.zeros(len(temps))
Cp_ave_O2_cal = np.zeros(len(temps))
Cp_ave_CO2_cal = np.zeros(len(temps))
Cp_ave_H2O_cal = np.zeros(len(temps))

for i in range(0,len(temps)):
    Cp_ave_Air_cal[i] = Cp_ave_Air_c(0,temps[i])
    Cp_ave_N2_cal[i] = Cp_ave_N2_c(0,temps[i])
    Cp_ave_O2_cal[i] = Cp_ave_O2_c(0,temps[i])
    Cp_ave_CO2_cal[i] = Cp_ave_CO2_c(0,temps[i])
    Cp_ave_H2O_cal[i] = Cp_ave_H2O_c(0,temps[i])

def Cp_ave_plot():
    plt.plot(temps,Cp_ave_Air(temps), label = 'Cp_ave_Air')
    plt.plot(temps,Cp_ave_N2s(temps), label = 'Cp_ave_N2*')
    plt.plot(temps,Cp_ave_N2(temps), label = 'Cp_ave_N2')
    plt.plot(temps,Cp_ave_O2(temps), label = 'Cp_ave_O2')
    plt.plot(temps,Cp_ave_CO2(temps), label = 'Cp_ave_CO2')
    plt.plot(temps,Cp_ave_H2O(temps), label = 'Cp_ave_H2O')

    #plt.plot(temps,Cp_ave_Air_cal, label = 'Cp_ave_Air_calc')
    plt.plot(temps,Cp_ave_N2_cal, label = 'Cp_ave_N2_calc')
    plt.plot(temps,Cp_ave_O2_cal, label = 'Cp_ave_O2_calc')
    plt.plot(temps,Cp_ave_CO2_cal, label = 'Cp_ave_CO2_calc')
    plt.plot(temps,Cp_ave_H2O_cal, label = 'Cp_ave_H2O_calc')
    
    plt.title('Cp_averages from 0 °C to t(°C)')
    plt.xlabel('Final Temperature [deg C]')
    plt.ylabel('Cp_ave in kJ/(kg K)')
    plt.legend(loc='best')
    plt.grid()
    plt.show()

