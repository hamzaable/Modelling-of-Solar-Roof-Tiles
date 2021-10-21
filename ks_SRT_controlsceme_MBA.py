# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:01:05 2021

@author: mariu
"""
import numpy as np
import scipy.interpolate
import pandas as pd

"______Import_Operating_strategies & Mass_Flow_Loss_values____________"

op_strategy = pd.read_excel(r'rauhigkeit.xlsx', sheet_name='Parameter')
op_strategy = op_strategy.drop(columns=['Einstrahlung [W/m2]', 'Umgebungstemperatur [C]', 'Windgeschwindigkeit [m/s]'])
op_strategy = op_strategy.assign(M = "")
op_strategy = op_strategy.rename(columns={'Unnamed: 0': 'Operating_Strategy', 'Volumemstrom [m3/h]': 'Volume_Flow_[m3/h]', 'M': 'Mass_Flow_[kg/s]'})

os_name = "5_1_dp{}"                                  # Creating variable for iterating through control stretegies 1 - 6
mass_flow_loss = pd.DataFrame({"SDP": ["SDP1", "SDP2", "SDP3", "SDP4",
                                           "SDP5", "SDP6", "SDP7", "SDP8",
                                           "SDP9", "SDP10", "SDP11", "SDP12"]})

for i in range(len(op_strategy)):
    m_flow_loss_temp = pd.read_excel(r'rauhigkeit.xlsx', sheet_name=os_name.format(str(i)), usecols = "A,I:J, L")
    op_strategy.loc[i, 'Mass_Flow_[kg/s]'] = m_flow_loss_temp.iloc[12]['m_dot']
    m_flow_loss_temp = m_flow_loss_temp.drop(range(12,28))
    mass_flow_loss.insert(i+1, str(os_name.format(str(i))),
                          m_flow_loss_temp["m_dot_leakage_h"] + m_flow_loss_temp["m_dot_leakage_v"])

mass_flow_loss = mass_flow_loss.select_dtypes(exclude='object').div(1.3).combine_first(mass_flow_loss)
first_c = mass_flow_loss.pop('SDP')                                                           # Reassemble Dataframe
mass_flow_loss.insert(0, 'SDP', first_c)

# mass_flow_loss = 0.001
# mass_flow_loss = None


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

i=0
mass_flow = 0.0320168





# iterate just one time + skip if mass flow is constant!! skip if mass flow is constant!


"""
function call: 
m_loss, ks_SRT, mass_flow_temp= interpolate_ks_mloss(i=i, 
                                                     op_strategy=op_strategy, 
                                                     os_name=os_name, 
                                                     mass_flow=mass_flow,
                                                     mass_flow_loss=mass_flow_loss
                                                     mass_flow_temp=mass_flow_temp):

function definition:
    def interpolate_ks_mloss(self,
                             i,
                             op_strategy,
                             os_name,
                             mass_flow,
                             mass_flow_temp):
    

"""
mflow = [0.064616108, 0.032016786, 0.02390436, 0.015828479, 0.009407934, 0.00306056]            # range of massflow values from operationg points DP3 to DP9 in kg/s
ks_mflow = [0.00018, 0.000225, 0.00024, 0.000295, 0.00035, 0.0021]                              # manually estimated ks values for Operating points DP3 to DP9 [-]

mass_flow_temp = None # from SDP

mflow = [0.064616108, 0.032016786, 0.02390436, 0.015828479, 0.009407934, 0.00306056]            # range of massflow values from operationg points DP3 to DP9 in kg/s
ks_mflow = [0.00018, 0.000225, 0.00024, 0.000295, 0.00035, 0.0021]                              # manually estimated ks values for Operating points DP3 to DP9 [-]

ks_interpolate = scipy.interpolate.interp1d(mflow, ks_mflow)
ks_interp = ks_interpolate(mass_flow) 
    
if (0.0646161 < mass_flow):
    raise ValueError('Mass Flow value exeeds maximum value of 200 m^3/h .')
elif (0.0646161 >= mass_flow >= 0.0320168):
    y2_str = os_name.format(str(8))
    y1_str = os_name.format(str(0))
    x2 = op_strategy[op_strategy['Operating_Strategy'].str.match(y2_str)].iloc[0][2]
    x1 = op_strategy[op_strategy['Operating_Strategy'].str.match(y1_str)].iloc[0][2]
elif (0.0320168 > mass_flow >= 0.0239044):
    y2_str = os_name.format(str(0))
    y1_str = os_name.format(str(5))
    x2 = op_strategy[op_strategy['Operating_Strategy'].str.match(y2_str)].iloc[0][2]
    x1 = op_strategy[op_strategy['Operating_Strategy'].str.match(y1_str)].iloc[0][2]
elif (0.0239044 > mass_flow >= 0.0158285):
    y2_str = os_name.format(str(5))
    y1_str = os_name.format(str(4))
    x2 = op_strategy[op_strategy['Operating_Strategy'].str.match(y2_str)].iloc[0][2]
    x1 = op_strategy[op_strategy['Operating_Strategy'].str.match(y1_str)].iloc[0][2]
elif (0.0158285 > mass_flow >= 0.00940793):
    y2_str = os_name.format(str(4))
    y1_str = os_name.format(str(9))
    x2 = op_strategy[op_strategy['Operating_Strategy'].str.match(y2_str)].iloc[0][2]
    x1 = op_strategy[op_strategy['Operating_Strategy'].str.match(y1_str)].iloc[0][2]
elif (0.00940793 > mass_flow >= 0.00621843):
    y2_str = os_name.format(str(9))
    y1_str = os_name.format(str(7))
    x2 = op_strategy[op_strategy['Operating_Strategy'].str.match(y2_str)].iloc[0][2]
    x1 = op_strategy[op_strategy['Operating_Strategy'].str.match(y1_str)].iloc[0][2]
elif (0.00621843 > mass_flow >= 0.00306056):
    y2_str = os_name.format(str(7))
    y1_str = os_name.format(str(6))
    x2 = op_strategy[op_strategy['Operating_Strategy'].str.match(y2_str)].iloc[0][2]
    x1 = op_strategy[op_strategy['Operating_Strategy'].str.match(y1_str)].iloc[0][2]
else:
    raise ValueError('Mass Flow value is below minimum value of 10 m^3/h .')
    
mass_flow_loss_interp = []
for k in range(len(mass_flow_loss)):
    y2 = mass_flow_loss[y2_str].iloc[k]
    y1 = mass_flow_loss[y1_str].iloc[k]
    y = round(y1 + ((y2-y1)/(x2-x1))*(mass_flow-x1), 9)     #later change in mass_flow[i]
    mass_flow_loss_interp.append(y)
    
mass_flow_loss_interp=pd.DataFrame(mass_flow_loss_interp)

mass_flow_temp = mass_flow

#return mass_flow_loss_interp, ks_interp, mass_flow_temp
