# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:55:55 2021

@author: marius bartkowski
"""
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

# Optical losses
tau = 0.93                                                                      # reflection loss
alpha = 0.90                                                                    # absorption loss
F_ta_en = tau * alpha                                                           # Null-Verlust-Effizienz F'(ta)en

# function definition
def func(x, c1, c2):
    return  F_ta_en * G - c1*x - c2*np.power(x, 2)                              # F_ta_en * G - c1*x - c2*x^2: See formula from Energiekonzept Documentation

# Input Data (later on with imported measurement values)
G = np.array([100, 250, 350, 500])                                              # Global Irradiance on collector surface (W/m^2)      
ta_2ms = np.array([8, 10, 12.5, 15])                                            # Ambient temperature at wind speed x (°C)
tm = np.array([10, 13.5, 15, 17.5])                                             # module temperature (°C)
dt = np.subtract(tm, ta_2ms)                                                    # delta t (°C)

# Data Points
yData = np.zeros(len(dt))                                                       # All zero because: P_out = 0 to simulate collector temperature in balance of losses and gains through solar irradiation

#Solve
InitialGuess = [1.0, 1.0]                                                       # initial guess necessary for the solver
coeff, pcov = curve_fit(func, dt, yData, InitialGuess)

print('Heat loss Coefficients:\n\nc1: ', round(coeff[0], 2), '(W/m*K)\nc2: ', round(coeff[1], 2), '(W/m*K^2)')

#Curve Fitting and Plotting
#xFit = np.arange(0, 30, 0.05)
#plt.plot(xFit, func(xFit, *coeff))

plt.plot(G, tm)
plt.xlabel('Einstrahlung (W/m^2)')
plt.ylabel('Modutemperatur (°C)')
plt.title('Measurement values of Module Temperature over different Irradiances')

