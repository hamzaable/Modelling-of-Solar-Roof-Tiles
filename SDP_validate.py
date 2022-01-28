# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 10:20:00 2022

@author: Marius Bartkowski
"""

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import datetime as dt
import pvlib

latitude = 50.934055
longitude = 6.990349
num_sdp_series = 96                                                             
num_sdp_parallel = 2

"_____________Data Imports_____________"
#Importing Measurement Data from Excel
os.chdir('C:/Users/mariu/Documents/GitHub/Modelling_of_Solar_Roof_Tiles/Modelling-of-Solar-Roof-Tile')
path = os.path.join("Imports", r'Datensatz_Komplett_2021-12-21.xlsx')

"__________Atmospheric conditions_____________"
print("\nImporting ...")
mdata_atmospheric_conditions = pd.read_excel(path, usecols = "A, H, J, L, AI, AN:AP", names=['Zeitstempel', 'wind_speed', 'temp', 'pr', 'T_inlet', 'DTI', 'DHI', 'GHI'])
mdata_atmospheric_conditions['GHI'][mdata_atmospheric_conditions['GHI'] <= 0] = 0
mdata_atmospheric_conditions['DHI'][mdata_atmospheric_conditions['GHI'] < mdata_atmospheric_conditions['DHI']] = mdata_atmospheric_conditions['GHI']
mdata_atmospheric_conditions['DNI'] = mdata_atmospheric_conditions["GHI"]-mdata_atmospheric_conditions["DHI"]



"__________Technical conditions_____________"
mdata_technical_conditions = pd.read_excel(path, usecols = "AF:AG, AM, CC, DG, EK, FN, FQ:RA, SA:SB")
Spannung_string_ohneMPPT = mdata_technical_conditions["Spannung_V_280"] + mdata_technical_conditions["Spannung_V_305"]+mdata_technical_conditions["Spannung_V_330"]+mdata_technical_conditions["Spannung_V_355"] 

sensor_name = "Spannung_V_{}" 
Spannung_string_mitMPPT = 0
sensor_count=0

for i in range(357, 463):    
    print(sensor_name.format(str(i)))
    
    try:
        Spannung_string_mitMPPT+= mdata_technical_conditions[sensor_name.format(str(i))]
        print(" voltage added")
        sensor_count+=1
    except:
        print(" passed\n")

print("valid voltage sensors: ", sensor_count) 

dcpower_string_ohneMPPT = Spannung_string_ohneMPPT * mdata_technical_conditions["Strom_A_475"]
dcpower_string_mitMPPT = Spannung_string_mitMPPT * mdata_technical_conditions["Strom_A_476"]

mdata =  mdata_atmospheric_conditions
mdata['th_heatflux_kW'] = mdata_technical_conditions.Th_Leistung_1_kW
mdata['mass_flow'] = mdata_technical_conditions.Massenstrom_1_kgh
mdata['mass_flow'] = mdata['mass_flow']/(3600*12)                               #conversion from kg/h to kg/s for one SRT string
mdata['HP_cons_elec_kW'] = mdata_technical_conditions.Strom_WP_  
mdata['Voltage_left_roof_part_V'] = Spannung_string_ohneMPPT
mdata['Current_left_roof_part_A'] = mdata_technical_conditions["Strom_A_475"]
mdata['DC-power_left_roof_part_W'] = dcpower_string_ohneMPPT
mdata['Voltage_Right_roof_part_V'] = Spannung_string_mitMPPT
mdata['Current_Right_roof_part_A'] = mdata_technical_conditions["Strom_A_476"]
mdata['DC-power_Right_roof_part_W'] = dcpower_string_mitMPPT

dcpower_string_ohneMPPT_sum = 0 
dcpower_string_mitMPPT_sum = 0
    
for i in range(len(mdata['DC-power_Right_roof_part_W'])-1):
    
    time_diff = mdata["Zeitstempel"][i+1] - mdata["Zeitstempel"][i]
         
    time_conversion_flux = time_diff.seconds * (60*60)  #1kW*Xsec = (1kW*Xsec*1min)/60sec = (1kW*1h*1min)/(60*60sec) = (1kW*1h)/(60*60) = 1kWh/(3.600)
    time_conversion_dc = time_diff.seconds/(60*60*1000) #1W*Xsec = (1W*Xsec*1min)/60sec = (1W*1h*1min)/(60*60sec) = (1W*1h)/(60*60) = 1kWh/(3.600.000)
    
    flux_sum = mdata['th_heatflux_kW'][i] * time_conversion_flux
    dcpower_string_ohneMPPT_sum += mdata['DC-power_left_roof_part_W'][i] * time_conversion_dc
    dcpower_string_mitMPPT_sum += mdata['DC-power_Right_roof_part_W'][i] * time_conversion_dc

print("DC-energy yield left roof part (without MPPT)", dcpower_string_ohneMPPT_sum, "kWh/Day")
print("DC-energy yield right roof part (with MPPT)", dcpower_string_mitMPPT_sum, "kWh/Day")
print("Heat Flux is: ", flux_sum, "kWh/Day")
print("normed Heat Flux is: ", flux_sum/(num_sdp_series*num_sdp_parallel*0.1), "kWh/m2*Day")           # 0.1 is module area

"______________Importing_AC-Power_Measurement_Data___________"

print("\nImporting ...")
path = os.path.join("Imports", r'ACPower_ohne_MPPT_21-12-21_15min.xlsx')
mdata_AC_power = pd.read_excel(path, usecols = "A, D", names=['Zeitstempel', 'AC-power_left_without_MPPT_kW'])
path = os.path.join("Imports", r'ACPower_MPPT_21-12-21_15min.xlsx')
mdata_AC_power['AC-power_right_with_MPPT_kW'] = pd.read_excel(path, usecols = "D")

time_conversion = 15/60 #1kW*15min = (1kW*15min*1h)/(60min)

print("AC-power without_MPPT: ", mdata_AC_power["AC-power_left_without_MPPT_kW"].sum()*time_conversion, " kWh")
print("\nAC-Power with MPPT:", mdata_AC_power["AC-power_right_with_MPPT_kW"].sum()*time_conversion, " kWh")

mdata.to_excel(os.path.join("Exports", r'MeasurementValues_calculated.xlsx'))


