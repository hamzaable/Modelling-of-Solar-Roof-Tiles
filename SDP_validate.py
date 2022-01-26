# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 10:20:00 2022

@author: Marius Bartkowski
- Class for importing Measurement Data 
"""

import os
import pandas as pd

class MeasurementDataImport():
    def ReturnMeasurementData(self):
    
    
        "_____________Data Imports_____________"
        #Importing Measurement Data from Excel
        #os.chdir('C:/Users/mariu/Documents/GitHub/Modelling_of_Solar_Roof_Tiles/Modelling-of-Solar-Roof-Tile')
        path = os.path.join("Imports", r'Datensatz_Komplett_2021-12-21.xlsx')
        
        "__________Atmospheric conditions_____________"
        print("\nImporting measurement Data for Atmospheric conditions ...")
        mdata_atmospheric_conditions = pd.read_excel(path, usecols = "A, H, J, L, AI, AN:AP", names=['Zeitstempel', 'wind_speed', 'temp', 'pr', 'T_inlet', 'DTI', 'DHI', 'GHI'])
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
        
        time_conversion = 10/(60*60*1000) #1W*10sec = (1W*10*1min)/60sec = (10W*1h*1min)/(60*60sec) = (10W*1h)/(60*60) = 10kWh/(3.600.000)
        
        print("DC-energy yield left roof part (without MPPT)", dcpower_string_ohneMPPT.sum()*time_conversion, "kWh")
        print("DC-energy yield right roof part (with MPPT)", dcpower_string_mitMPPT.sum()*time_conversion, "kWh")
        
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
        
        "______________Importing_AC-Power_Measurement_Data___________"
        
        print("\nImporting measurement Data  for AC Output...")
        path = os.path.join("Imports", r'ACPower_ohne_MPPT_21-12-21_15min.xlsx')
        mdata_AC_power = pd.read_excel(path, usecols = "A, D", names=['Zeitstempel', 'AC-power_left_without_MPPT_kW'])
        path = os.path.join("Imports", r'ACPower_MPPT_21-12-21_15min.xlsx')
        mdata_AC_power['AC-power_right_with_MPPT_kW'] = pd.read_excel(path, usecols = "D")
        
        time_conversion = 15/60 #1kW*15min = (1kW*15min*1h)/(60min)
        
        print("\nAC-Power with MPPT:", mdata_AC_power["AC-power_right_with_MPPT_kW"].sum()*time_conversion, " kWh")
        print("AC-power without_MPPT: ", mdata_AC_power["AC-power_left_without_MPPT_kW"].sum()*time_conversion, " kWh")
        
        return mdata_AC_power, mdata
        
        

