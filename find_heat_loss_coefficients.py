# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:55:55 2021

@author: marius bartkowski
"""
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
import seaborn as sns
import numpy as np
import pandas as pd
import os
import math
from  matplotlib.ticker import FuncFormatter

# Path handling
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))
#Change Path of directory manually
project_path = 'C:/Users/mariu/Documents/GitHub/Modelling_of_Solar_Roof_Tiles/Modelling-of-Solar-Roof-Tile'
os.chdir(project_path)
# Change the current working directory automaticly
#os.chdir(os.path.abspath('find_heat_loss_coefficients.py'))
print("\nChanged working directory to:\n {0}".format(os.getcwd()))


# Import Data
print("\nImporting ...\n\n")
can_data = pd.read_excel(os.path.join("Imports", r'can_data.xlsx'), usecols = "BI")


weather_data = pd.read_excel(os.path.join("Imports", r'weather_data.xlsx'), usecols = "A:C, E, I")
meta = pd.concat([weather_data, can_data], axis=1) 
meta = weather_data.reset_index().merge(can_data.reset_index()).drop(columns='index').rename(columns={'Temperature °C': 'Temperature ambient °C',
                                                                                                      'Temperatur_1_C_53': 'Temperature 1 module C53',
                                                                                                      'Global Irradiance W/m2': 'Global_Irradiance_W_m2'
                                                                                                      })

meta = meta[meta.Global_Irradiance_W_m2 > 30]                                   #filter out every datapoint below Irradiance of 30 W/m2 due to unreasonable values

"______________Datetime calculations_____________________"

# Really ugly made data selection
meta.reset_index(inplace=True, drop=True) 
meta_11_10 = meta[0:2624]                                                       #select data for 11.11.21 between 9:26:20 and 16:39:40
meta_11_11 = meta[3056:5629]                                                    #select data for 11.11.21 between 9:26:20 and 16:39:40

meta_ws2 = meta_11_10.append(meta_11_11)                                        #merge both dataframes 

meta_ws2 = meta_ws2[(meta_ws2["Windspeed m/s"] < 2.05 ) & (meta_ws2["Windspeed m/s"] > 1.95 )]    #select data with Windspeed 2 m/s +- 0,5 m/s - Originally 

# Optical losses
tau = 0.93                                                                      # reflection loss
alpha = 0.90                                                                    # absorption loss
F_ta_en = tau * alpha                                                           # Null-Verlust-Effizienz F'(ta)en

# function definition
def func(x, c1, c2):
    return  F_ta_en * G - c1*x - c2*np.power(x, 2)                              # F_ta_en * G - c1*x - c2*x^2: See formula from Energiekonzept Documentation

# Input Data (imported measurement data)
G = np.array(meta_ws2["Global_Irradiance_W_m2"])                                # Global Irradiance on collector surface (W/m^2)      
ta_2ms = np.array(meta_ws2["Temperature ambient °C"])                           # Ambient temperature at wind speed x (°C)
tm = np.array(meta_ws2["Temperature 1 module C53"])                             # module temperature (°C)
dt = np.subtract(tm, ta_2ms)                                                    # delta t (°C)
ws = np.array(meta_ws2["Windspeed m/s"])                                        # Windspeed (m/s)                

# Data Points
yData = np.zeros(len(dt))                                                       # All zero because: P_out = 0 to simulate collector temperature in balance of losses and gains through solar irradiation

# Solve
InitialGuess = [1.0, 1.0]                                                       # initial guess necessary for the solver
coeff, pcov = curve_fit(func, dt, yData, InitialGuess)

print('Heat loss Coefficients via measurement data:\n\nc1: ', round(coeff[0], 2), '(W/m*K)\nc2: ', round(coeff[1], 2), '(W/m*K^2)')


plt.figure()
ax = sns.regplot(x=G, y=tm, scatter_kws={"color": "black"}, line_kws={"color": "red"})
ax.set_yticklabels(ax.get_yticks(), size = 14)
ax.set_xticklabels(ax.get_xticks(), size = 14)
plt.xlabel('total Irradiance on the collector plane (W/m2)', fontsize=16)
plt.ylabel('Module Temperature (°C)', fontsize=16)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: int(y)))
#plt.title('Measurement values of Module Temperature over different Irradiances')


# Calculating module temperature + c1 & c2 via Christian Brosig's approach

a = -2.98                                                                       # empirical coefficient : upper limit of module temperature at high irradiances
b = -0.0471                                                                     # empirical coefficient : rate of sinking moduel tempeature at rising windspeed

tm_sandia = np.empty(len(ta_2ms))

for i in range(len(ta_2ms)) :
    tm_sandia[i] = G[i] * (math.exp(a + (b * ws[i]))) + ta_2ms[i]               # source: Documentation energiekonzept
    
dt_sandia = np.subtract(tm_sandia, ta_2ms)

plt.figure()
ax = sns.regplot(x=G, y=tm_sandia)
plt.xlabel('Irradiance (W/m2)')
plt.ylabel('Module Temperature (°C)')
plt.title('Calculated values for Module Temperature over different Irradiances via Sandia model')

# Solve
InitialGuess = [1.0, 1.0]                                                       # initial guess necessary for the solver
coeff_tm_Sandia, pcov = curve_fit(func, dt_sandia, yData, InitialGuess)

print('\n\nHeat loss Coefficients via module temperature by Sandia model:\n\nc1: ', round(coeff_tm_Sandia[0], 2), '(W/m*K)\nc2: ', round(coeff_tm_Sandia[1], 2), '(W/m*K^2)')

"___________Collector efficiency____________________"

col_name = "η SRT ({} W/m2)"

G_ref = 1000 #in [W/m^2]
dt_eff = np.arange(0, 161, 1)

df_eta_SRT = pd.DataFrame(data= dt_eff, columns={'dt [°C]'})  

while(G_ref>0):
    list_eta_SRT  = []
    for i in range(len(dt_eff)):
        eta_SRT = F_ta_en - (round(coeff[0], 2) * (dt_eff[i]/G_ref)) - (round(coeff[1], 2)*(np.power((dt_eff[i]), 2)/G_ref))
        list_eta_SRT.append(eta_SRT)
        
    df_eta_SRT[col_name.format(str(G_ref))] = list_eta_SRT
    G_ref-=200

plt.figure(figsize=(8, 6))
cmap = cm.get_cmap('viridis')
ax = df_eta_SRT.plot(x="dt [°C]", y=["η SRT (1000 W/m2)", "η SRT (800 W/m2)", "η SRT (600 W/m2)", "η SRT (400 W/m2)", "η SRT (200 W/m2)"],linewidth=3.0, cmap=cmap)
params = {'mathtext.default': 'regular' } 
plt.rcParams.update(params)
plt.xlabel('$T_{m} - T_{a}  [°C]$', fontsize=15)
plt.ylabel('collector efficiency η [-]', fontsize=15)
plt.legend(fontsize=12)
ax.set_xlim(0,60)
ax.set_yticklabels(ax.get_yticks(), size = 14)
ax.set_xticklabels(ax.get_xticks(), size = 14)
ax.set_ylim(0,1.0)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: round(float(y), 1)))




