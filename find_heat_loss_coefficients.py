# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:55:55 2021

@author: marius bartkowski
"""
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import os
#import time
import datetime as dt
import math

# Path handling
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))
project_path = 'C:/Users/mariu/Documents/GitHub/Modelling_of_Solar_Roof_Tiles/Modelling-of-Solar-Roof-Tile'
# Change the current working directory
os.chdir(project_path)
print("\nChanged working directory to:\n {0}".format(os.getcwd()))


# Import Data
print("\nImporting ...\n\n")
can_data = pd.read_excel(r'can_data.xlsx', usecols = "BI")         

weather_data = pd.read_excel(r'weather_data.xlsx', usecols = "A:C, E, I") 
meta = pd.concat([weather_data, can_data], axis=1) 
meta = weather_data.reset_index().merge(can_data.reset_index()).drop(columns='index').rename(columns={'Temperature °C': 'Temperature ambient °C',
                                                                                                      'Temperatur_1_C_53': 'Temperature 1 module C53',
                                                                                                      'Global Irradiance W/m2': 'Global_Irradiance_W_m2'
                                                                                                      })

meta = meta[meta.Global_Irradiance_W_m2 > 30]                                   #filter out every datapoint below Irradiance of30 W/m2 due to unreasonable values

"______________Datetime calculations_____________________"
#Choosed an easier way because its not working
"""
nineteen_seventy = time.strptime('01-01-70', '%d-%m-%y')

print(meta_test['Time']) #meta_test.Global_Irradiance_W_m2 > 30])
print(meta_test[meta_test.Time])

con_time = pd.to_datetime(meta_test['Time'])
con_time = con_time.dt.time

con_time.dtypes

time_str = str(meta_test['Time'])
meta_test.dtypes
"""

# Really ugly made data selection
meta.reset_index(inplace=True, drop=True) 
meta_11_10 = meta[0:2624]                                                       #select data for 11.11.21 between 9:26:20 and 16:39:40
meta_11_11 = meta[3056:5629]                                                    #select data for 11.11.21 between 9:26:20 and 16:39:40

meta_ws2 = meta_11_10.append(meta_11_11)                                        #merge both dataframes 

meta_ws2 = meta_ws2[(meta_ws2["Windspeed m/s"] < 2.5 ) & (meta_ws2["Windspeed m/s"] > 1.5 )]    #select data with Windspeed 2 m/s +- 0,5 m/s

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

# Data Points
yData = np.zeros(len(dt))                                                       # All zero because: P_out = 0 to simulate collector temperature in balance of losses and gains through solar irradiation

#Solve
InitialGuess = [1.0, 1.0]                                                       # initial guess necessary for the solver
coeff, pcov = curve_fit(func, dt, yData, InitialGuess)

print('Heat loss Coefficients:\n\nc1: ', round(coeff[0], 2), '(W/m*K)\nc2: ', round(coeff[1], 2), '(W/m*K^2)')

#Curve Fitting and Plotting - Not working yet
"""
n = len(dt)
y = np.empty(len(dt)) 
for i in range(n) :
    y[i] = func(G[i], coeff[0], coeff[1])
    
# function set to tm = ...
# y[i] = (math.sqrt((coeff[0]^2)+4*coeff[1]*F_ta_en*G[i])-coeff[0]+(2*coeff[1]*ta_2ms[i]))/(2*coeff[1])
"""

plt.plot(G, tm, 'k.')
plt.xlabel('Einstrahlung (W/m^2)')
plt.ylabel('Modutemperatur (°C)')
plt.title('Measurement values of Module Temperature over different Irradiances')



