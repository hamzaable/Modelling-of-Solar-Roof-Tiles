# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:01:05 2021

@author: mariu
"""

"________________ks_SRT in dependency of weather data_________________________________"
ks_SRT = 0.00007 # one year calculation overall ks value for pr 250 Pa



ghi = 805
temp_amb = 20
wind_amb = 735

if (ghi > 800) and (temp_amb <= 25) and (wind_amb <= 1):
    ks_SRT_in = 0
elif (ghi > 800) and (temp_amb > 25) and (wind_amb <= 1):
    ks_SRT_in = 1
elif (ghi <= 800) and (temp_amb <= 25) and (wind_amb <= 1):
    ks_SRT_in = 2
elif (ghi > 800) and (temp_amb <= 25) and (wind_amb > 1):
    ks_SRT_in = 3
else:
    ks_SRT_in = ks_SRT

print("Irradiance: ", ghi, "Temp: ", temp_amb, "Wind ", wind_amb, "ks --> ", ks_SRT_in) 



"________________mass_flow_loss in dependency of mass Flow value__________________________"
ks_SRT = 0.00007 # one year calculation overall ks value for pr 250 Pa



ghi = 805
temp_amb = 20
wind_amb = 735

if (ghi > 800) and (temp_amb <= 25) and (wind_amb <= 1):
    ks_SRT_in = 0
elif (ghi > 800) and (temp_amb > 25) and (wind_amb <= 1):
    ks_SRT_in = 1
elif (ghi <= 800) and (temp_amb <= 25) and (wind_amb <= 1):
    ks_SRT_in = 2
elif (ghi > 800) and (temp_amb <= 25) and (wind_amb > 1):
    ks_SRT_in = 3
else:
    ks_SRT_in = ks_SRT

print("Irradiance: ", ghi, "Temp: ", temp_amb, "Wind ", wind_amb, "ks --> ", ks_SRT_in) 

