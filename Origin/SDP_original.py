# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 22:15:04 2020

@author: Vaishnavi Phadke
"""
import pandas as pd
import numpy as np
import os

##########
from PVLIB_model_original import Photovoltaic
from TESPy_model_original import SDP_sucking
##########
start = "01-01-{} 00:00".format(str(2019))
end = "31-12-{} 23:00".format(str(2019))
naive_times = pd.date_range(start= start, end=end, freq='1h')

latitude = 50.9375
longitude = 6.9603
name='Cologne'
altitude=  121
timezone= 'Etc/GMT+2'

"_____________Data Imports_____________"
"Hourly Weather Data (DNI , GHI , DHI , temp_air , wind speed and pressure)"
wd_path = "C:/Users/mariu/Documents/GitHub/Modelling_of_Solar_Roof_Tiles/Modelling-of-Solar-Roof-Tile/Origin"

dwd_data = pd.read_excel(os.path.join(wd_path, r'Timeseries_JRC_PVGIS_TH_Koeln.xlsx'))  # Hourly Weather Data (DNI , GHI , DHI , temp_air , wind speed and pressure)

pv_data = pd.DataFrame(index=dwd_data.index, columns=["dni","ghi",
                                                     "dhi",
                                                     "temp_air",
                                                     "wind_speed", "pressure"])

pv_data["DateTimeIndex"] = dwd_data.date
pv_data["DateTimeIndex"] = pd.to_datetime(pv_data["DateTimeIndex"])
pv_data["dni"] = dwd_data.irradiance_dir
pv_data["dhi"] = dwd_data.irradiance_diff
pv_data["ghi"] = dwd_data.poa_global
pv_data["temp_air"] = dwd_data.temp
pv_data["wind_speed"] = dwd_data.wind_speed
pv_data["pressure"] = dwd_data.pr

house_data_read = pd.read_excel(os.path.join(wd_path, r'house_demand.xlsx'))
house_data = pd.DataFrame(index=house_data_read.index, columns=["elec_cons","thermal_cons"])

house_data["DateTimeIndex"] = house_data_read.date
house_data["DateTimeIndex"] = pd.to_datetime(house_data["DateTimeIndex"])
house_data["elec_cons"] = house_data_read.elec_cons
house_data["thermal_cons"] = house_data_read.thermal_cons


num_sdp_series=12
num_sdp_parallel=16     # Only 12 for tespy calculations
module_number = num_sdp_series*num_sdp_parallel

"""
#####
#Thermal initialization
#####
sdp = SDP_sucking(sdp_in_parallel=num_sdp_parallel,
                  sdp_in_series=num_sdp_series)


sdp.init_sdp(ambient_temp=-4,
             absorption_incl=300,
             inlet_temp=-4,
             mass_flow=0.0320168,
             # zeta=4e6,
             m_loss=0.001,
             print_res=False)

######
##Electrical Yeild
######

"""
df = []
df1= []

for i in pv_data.index[0:8759]:
    time = pv_data.DateTimeIndex[i]
    temp_amb = pv_data.temp_air[i]
    wind_amb = pv_data.wind_speed[i]
    ghi = pv_data.ghi[i]
    dni = pv_data.dni[i]
    dhi = pv_data.dhi[i]
    
    electrical_yield= Photovoltaic(latitude=latitude, longitude=longitude, altitude=altitude, timezone=timezone, time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i], dhi=pv_data.dhi[i], temp_amb=pv_data.temp_air[i], wind_amb=pv_data.wind_speed[i],pressure=pv_data.pressure[i])
    df1 = [i, time, temp_amb, round(electrical_yield.annual_energy*module_number,2), int(electrical_yield.effective_irradiance), round(electrical_yield.tcell.item(), 2)]
    df.append(df1)

column_values = ["Index","Time","Tamb","Power", "Effective_Irradiance", "TModule"]
electrical_data = pd.DataFrame(data=df, columns = column_values)
electrical_data.fillna(0)
electrical_data.loc['Total'] = electrical_data.select_dtypes(np.number).sum()
pd.set_option('display.max_colwidth', 40)
print(electrical_data)


electrical_data.to_excel(r'Results1.xlsx') 

#####
#thermal Yeild
#####
"""
y=0
df = []
df1= []
for i in pv_data.index[0:8760]:
    time = pv_data.DateTimeIndex[i]    
    E_sdp = (0.93*(electrical_data.Effective_Irradiance[i]*0.10)-(electrical_data.Power[i]))/0.10 
    Tamb = pv_data.temp_air[i]
        
    if (E_sdp==0):
            #in deg celcius
            t_out_init = Tamb

            #Watt
            p_fan_init = 0
 
            #kg/sec
            m_out_init = 0 
            
    else:
            t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(ambient_temp=pv_data.temp_air[i], absorption_incl=E_sdp, inlet_temp=pv_data.temp_air[i], mass_flow=1, print_res=False)
                        
    t_out =  t_out_init
    p_fan =  p_fan_init
    m_out =  m_out_init
    flux = round((m_out*(t_out-Tamb)/(num_sdp_series*num_sdp_parallel*0.10)),2)
    
    elec_parameter = (house_data.elec_cons[i]+p_fan)<(electrical_data.Power[i]*num_sdp_parallel*num_sdp_series) 
    thermal_parameter = (house_data.thermal_cons[i]<flux)
    
    if(E_sdp==0):
        status = "System Off"
    
    elif (elec_parameter==True)and(thermal_parameter==True):
        status = "Good"
        
    elif (elec_parameter==True)and(thermal_parameter==False):
        status = "Increase Thermal Production (increase fan_power)"
        # t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(ambient_temp=pv_data.temp_air[i], absorption_incl=pv_data.ghi[i], inlet_temp=pv_data.temp_air[i], mass_flow=0.42, print_res=False)
        # t_out =  t_out_init
        # p_fan =  p_fan_init
        # m_out =  m_out_init
        # flux = round((m_out*(t_out-Tamb)/(num_sdp_series*num_sdp_parallel*0.10)),2)
        
        # elec_parameter = (house_data.elec_cons[i]+p_fan)<(electrical_data.Power[i]*num_sdp_parallel*num_sdp_series) 
        # thermal_parameter = (house_data.thermal_cons[i]<flux)
        
        # elec_parameter = (house_data.elec_cons[i]+p_fan)<(electrical_data.Power[i]*num_sdp_parallel*num_sdp_series) 
        # thermal_parameter = (house_data.thermal_cons[i]<flux)
        
    elif (elec_parameter==False)and(thermal_parameter==True):
        status = "Decrease Thermal Production (reduce fan_power)" 
        # t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(ambient_temp=pv_data.temp_air[i], absorption_incl=pv_data.ghi[i], inlet_temp=pv_data.temp_air[i], mass_flow=0.42, print_res=False)
        # t_out =  t_out_init
        # p_fan =  p_fan_init
        # m_out =  m_out_init
        # flux = round((m_out*(t_out-Tamb)/(num_sdp_series*num_sdp_parallel*0.10)),2)
        
        # elec_parameter = (house_data.elec_cons[i]+p_fan)<(electrical_data.Power[i]*num_sdp_parallel*num_sdp_series) 
        # thermal_parameter = (house_data.thermal_cons[i]<flux)
        
        # elec_parameter = (house_data.elec_cons[i]+p_fan)<(electrical_data.Power[i]*num_sdp_parallel*num_sdp_series) 
        # thermal_parameter = (house_data.thermal_cons[i]<flux)
        
    else:   
        status = "System Off"
        
    df1 = [i, time, Tamb, round(E_sdp ,2), t_out, p_fan, m_out, flux, status, elec_parameter, thermal_parameter]
    df.append(df1)
print(df)

column_values = ["Index","Time", "Tamb", "E_sdp_eff","T_out", "P_fan", "M_out", "HeatFlux", "status", "Elec_demand_met", "Heat_demand_met"]
thermal_data = pd.DataFrame(data=df, columns = column_values)
thermal_data.loc['Total'] = thermal_data.select_dtypes(np.number).sum()
pd.set_option('display.max_colwidth', 8)
print(thermal_data)

#Efficiency wrt effective E sdp for solar
Efficiency = thermal_data.loc["Total","HeatFlux"]/thermal_data.loc["Total","E_sdp_eff"]
print("Efficiency wrt Effective Irradiance:" , round(Efficiency*100,2),"%")


complete_data=pd.merge(electrical_data,thermal_data)
complete_data.to_excel(r'Result2.xlsx')
"""
