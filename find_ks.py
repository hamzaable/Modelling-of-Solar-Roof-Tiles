# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 17:12:54 2021

@author: mariu
"""

# -*- coding: utf-8 -*-
#  Test cooling effect
import pandas as pd
import time
import datetime
from tespy.components import (Sink, Source, Valve, Merge,
                              Subsystem, Compressor,
                              SolarCollector)

##########
#from PVLIB_model import Photovoltaic
#from PVLIB_model import cellTemperature
from TESPy_model import SDP_sucking
#from TESPy_model import interpolate_ks_mloss

##########
start = "01-01-{} 00:00".format(str(2019))
end = "31-12-{} 23:00".format(str(2019))
naive_times = pd.date_range(start=start, end=end, freq='1h')
"________Commands & Control______________________"
find_ks_value = 0                                                               # If 0: ks values are not iterativly found, If 1: ks value will be found for each operating point

"________Location Parameters___________"
latitude = 50.9375
longitude = 6.9603
name = 'Cologne'
altitude = 121
timezone = 'Etc/GMT+2'

"______Photovoltaic Parameters_________"

albedo = 0.20                                                                   # Ground reflection albedo factor 0.20 (Beton) --> K. Mertens- Photovoltaik S.52
a_r = 0.14                                                                      # Spectral Corrections factor for different module glasses
irrad_model = 'haydavies'                                                       # Model for Irradiation calculation. Choose from: 'isotropic', 'klucher', 'haydavies', 'reindl', 'king', 'perez'
m_azimut = 180                                                                  # Module Azimut (Ausrichtung) [°]
m_tilt = 39                                                                     # Module tilt (Neigung) [°]

"""
        The module dict defined below is important for effective irradiance 
        calculations as well as the fitting of the IV-Curve via the fit-CEC-SAM 
        function, which results are used to calculate the PV-power output via 
        the single-diod-model. Within the module dict follwing specs can be 
        defined by the user:
        Specs:
        -----------
            A0 : float
                Airmass coefficient a0
            A1 : float
                Airmass coefficient a1
            A2 : float
                Airmass coefficient a2
            A3 : float
                Airmass coefficient a3
            A4 : float
                Airmass coefficient a4
            Aisc : float
                Temperature coefficient of short circuit current [A/C]
            Area : float
                Area of one cell [sqm]

            Bvoco : float
                Temperature coefficient of open circuit voltage [V/C]

            Cells_in_Series : int
                Number of cells in series [-]

            Impo : float
                Current at maximum power point [A]

            Isco : float
                Short circuit current [A]

            Material : str
                the length (or height) of the sdp in m

            Vintage : int
                Year of Construction

            Vmpo : float
                Voltage at maximum power point [V]

            Voco : float
                Open circuit voltage [V]

            celltype : str
                Value is one of ‘monoSi’, ‘multiSi’, ‘polySi’, ‘cis’, ‘cigs’, ‘cdte’, ‘amorphous’

            gamma_pmp : float
                Temperature coefficient of power at maximum point point [%/C]

"""
# ======Module Parameters=======================================================
# #ar=0.14
module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8,
          "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, "Aisc": 0.0010, "Bvoco": -0.0158,
          "gamma_pmp": -0.3792, "A0": 0.9645, "A1": 0.02753, "A2": -0.002848, "A3": -0.0001439,
          "A4": 0.00002219}
# =============================================================================

"_____________Data Imports_____________"
"Hourly Weather Data (DNI , GHI , DHI , temp_air , wind speed and pressure)"

dwd_data = pd.read_excel(r'Timeseries_JRC_PVGIS_TH_Koeln.xlsx')  # Hourly Weather Data (DNI , GHI , DHI , temp_air , wind speed and pressure)

pv_data = pd.DataFrame(index=dwd_data.index, columns=["dni", "ghi",
                                                      "dhi",
                                                      "temp_air",
                                                      "wind_speed", "pressure"])
"__________Assign Headers to pv_data Dataframe________________"

pv_data["DateTimeIndex"] = dwd_data.date
pv_data["DateTimeIndex"] = pd.to_datetime(pv_data["DateTimeIndex"])
pv_data["dni"] = dwd_data.irradiance_dir                                        # direct Solar irradition (horizontal)
pv_data["dhi"] = dwd_data.irradiance_diff                                       # diffuse solar radiation (horizontal)
pv_data["ghi"] = dwd_data.poa_global                                            # global irrradiation     (horizontal)
pv_data["temp_air"] = dwd_data.temp
pv_data["wind_speed"] = dwd_data.wind_speed
pv_data["pressure"] = dwd_data.pr

"________Hourly house demand (elec_cons , thermal_cons)________"
house_data_read = pd.read_excel(r'house_demand.xlsx')
house_data = pd.DataFrame(index=house_data_read.index, columns=["elec_cons", "thermal_cons"])

"_____________Assign Headers to pv_data Dataframe______________"
house_data["DateTimeIndex"] = house_data_read.date
house_data["DateTimeIndex"] = pd.to_datetime(house_data["DateTimeIndex"])
house_data["elec_cons"] = house_data_read.elec_cons
house_data["thermal_cons"] = house_data_read.thermal_cons

"______MassFlow Loss Import_________"
"""
# Import Mass flow loss table
mass_flow_loss = pd.read_excel(r'CFD_Daten_July.xlsx', sheet_name=1)  # mass flow losses for each SRT string in [kg/s]
mass_flow_loss = mass_flow_loss.drop(['position', 'Massenstrom_[kg/s]'], axis=1)  # Erase unnecessary columns
mass_flow_loss = mass_flow_loss.dropna()  # Erase Row with Nan values
mass_flow_loss = mass_flow_loss.reset_index(drop=True)  # Reset index
mass_flow_loss = mass_flow_loss.select_dtypes(exclude='object').div(1.3).combine_first(
    mass_flow_loss)  # Divide through factor 1.3
first_c = mass_flow_loss.pop('undichtigkeit')  # Reassemble Dataframe
mass_flow_loss.insert(0, 'undichtigkeit', first_c)
"""

"______Import_Operating_strategies & Mass_Flow_Loss_values____________"

#Operating Strategies
op_strategy = pd.read_excel(r'rauhigkeit.xlsx', sheet_name='Parameter')
op_strategy = op_strategy.assign(M = "", L = "")
op_strategy = op_strategy.rename(columns={'Unnamed: 0': 'Operating_Strategy',
                                          'Einstrahlung [W/m2]': 'Irradiance_[W/m2]',
                                          'Umgebungstemperatur [C]': 'Ambient_temperature [°C]',
                                          'Windgeschwindigkeit [m/s]': 'wind_speed [m/s]',
                                          'Volumemstrom [m3/h]': 'Volume_Flow_[m3/h]', 
                                          'M': 'Mass_Flow_[kg/s]',
                                          'L': 'Cumulative_Pressure_Drop_[Pa]'
                                          })
op_strategy = op_strategy.reindex(columns=['Operating_Strategy','Volume_Flow_[m3/h]','Mass_Flow_[kg/s]','Cumulative_Pressure_Drop_[Pa]','Irradiance_[W/m2]','Ambient_temperature [°C]','wind_speed [m/s]'])

os_name = "5_1_dp{}" # Creating variable for iterating through control stretegies 

# Mass Flow leakage
mass_flow_loss = pd.DataFrame({"SDP": ["SDP1", "SDP2", "SDP3", "SDP4",
                                           "SDP5", "SDP6", "SDP7", "SDP8",
                                           "SDP9", "SDP10", "SDP11", "SDP12"]})

for i in range(len(op_strategy)):
    m_flow_loss_temp = pd.read_excel(r'rauhigkeit.xlsx', sheet_name=os_name.format(str(i)), usecols = "A,I:J, L, R")
    op_strategy.loc[i, 'Mass_Flow_[kg/s]'] = m_flow_loss_temp.iloc[12]['m_dot']
    op_strategy.loc[i, 'Cumulative_Pressure_Drop_[Pa]'] = m_flow_loss_temp.iloc[27]['V_dot_leakage_h']
    m_flow_loss_temp = m_flow_loss_temp.drop(range(12,28))
    mass_flow_loss.insert(i+1, str(os_name.format(str(i))),
                          m_flow_loss_temp["m_dot_leakage_h"] + m_flow_loss_temp["m_dot_leakage_v"])

mass_flow_loss = mass_flow_loss.select_dtypes(exclude='object').div(1.3).combine_first(mass_flow_loss)      # all values for the leakge mass flow in the valves has to be divided by 1.3 
first_c = mass_flow_loss.pop('SDP')                                                                         # Reassemble Dataframe
mass_flow_loss.insert(0, 'SDP', first_c)

# mass_flow_loss = 0.001                                                        # assign one unitary value for the mass flow leakage in each valve
# mass_flow_loss = None

"______TESPy Model Parameters_________"

num_sdp_series = 12                                                             # Changed from 12 to 2 for test purpose
num_sdp_parallel = 12                                                           # Changed from 38 to 1 for test purpose
                                                                                # ks/roughness value for one SRT, used in design mode to calculate the pressure drop. ks_SRT values for off design mode are calculated
p_amb=1.01325                                                                   # Atmospheric pressure [Bar]
mass_flow = 0.0320168                                                           # Can be one value or string (from measurement data later on). IMPORTANT: This mass flow value applies for one String of 12 SRTs and is not the mass flow delivered by the fan for the whole SRT plant!
                                                             
# Allowed Value range for mass flow is:
# 0.0646 to 0.00306 kg/s (due to interpolation boundaries)

mass_flow_temp = None                                                           #Can be ignored

#####
# Thermal initialization
#####
sdp = SDP_sucking(sdp_in_parallel=num_sdp_parallel,
                  sdp_in_series=num_sdp_series,
                  #ks=ks
                  )

sdp.init_sdp(ambient_temp=-4,
             absorption_incl=300,
             inlet_temp=-4,
             mass_flow=op_strategy.iloc[0][2],
             #zeta=4e6,
             m_loss=mass_flow_loss,
             print_res=False)

"____Find_ks_values______"

t_start = time.time()
ks_result = []
iter_global=0

for i in range(len(op_strategy)): 
    iter_t = 0
    
    if op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][2] > 0.015:
        ks_SRT = 0
        p_deviation = -1
        while p_deviation <= -0.2:
            iter_t += 1
            ks_SRT += 0.000005                                                      #0.000005 --> originally         
            # for i in tqdm(pv_data.index[3000:3120]):      
            t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(
                      ambient_temp=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][5],
                      absorption_incl=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][4],
                      inlet_temp=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][5],
                      mass_flow=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][2],
                      print_res=False,
                      ks_SRT=ks_SRT,
                      m_loss_offdesign=pd.DataFrame(mass_flow_loss[str(os_name.format(str(i)))]).set_axis(['5_1_dpx'], axis=1, inplace=False)
                      )  
        
            p_Valve = []
            
            for comp in sdp.nw.comps['object']:
                if isinstance(comp, Valve):
                    p_Valve.append(comp.pr.val)
                    
            p_diff = round((p_Valve[11] - p_Valve[0])*100000, 2)
            p_deviation = round(op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][3]-p_diff, 2)
            print("ks value: ", round(ks_SRT, 7), "\nCalculated Pressure difference is: ", p_diff, "Pa\nDeviation to CFD values: ", p_deviation, "Pa\n")       
    
    else:
        p_deviation = 1
        ks_SRT = 0.0009
        while p_deviation >= 0.2:
            iter_t += 1
            ks_SRT -= 0.000005         
            # for i in tqdm(pv_data.index[3000:3120]):      
            t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(
                      ambient_temp=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][5],
                      absorption_incl=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][4],
                      inlet_temp=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][5],
                      mass_flow=op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][2],
                      print_res=False,
                      ks_SRT=ks_SRT,
                      m_loss_offdesign=pd.DataFrame(mass_flow_loss[str(os_name.format(str(i)))]).set_axis(['5_1_dpx'], axis=1, inplace=False)
                      )  
        
            p_Valve = []
            
            for comp in sdp.nw.comps['object']:
                if isinstance(comp, Valve):
                    p_Valve.append(comp.pr.val)
                    
            p_diff = round((p_Valve[11] - p_Valve[0])*100000, 2)
            p_deviation = round(op_strategy[op_strategy['Operating_Strategy'].str.match(os_name.format(str(i)))].iloc[0][3]-p_diff, 2)
            print("ks value: ", round(ks_SRT, 7), "\nCalculated Pressure difference is: ", p_diff, "Pa\nDeviation to CFD values: ", p_deviation, "Pa\n")
        
    ks_result.append(                                                                                         
                    {
                        'Operating Strategy': os_name.format(str(i)),                                                               
                        'ks [-]': round(ks_SRT, 7),                                                                      
                        'At Deviation to real value [Pa]': round(p_deviation, 2) 
                    }
                    )
    
    print("\n\nks value for strategy", os_name.format(str(i)), "is ", round(ks_SRT, 7), "\nContinue...\n\n")
    iter_global += iter_t
    
t_end = time.time()
print("Calculation time:", str(datetime.timedelta(seconds=t_end - t_start)), "[hh:mm:ss]\nIterations: ", iter_global, "\n\n")
ks_result = pd.DataFrame(ks_result)
print(ks_result, "\nResults stored in ks_result.")

