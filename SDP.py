# -*- coding: utf-8 -*-
#  Test cooling effect
import os

import pandas as pd
import numpy as np

from tqdm import tqdm

##########
from PVLIB_model import Photovoltaic
from PVLIB_model import cellTemperature
from TESPy_model import SDP_sucking
from heat_pump import HeatPump
#from TESPy_model import interpolate_ks_mloss

##########
start = "01-01-{} 00:00".format(str(2019))
end = "31-12-{} 23:00".format(str(2019))
naive_times = pd.date_range(start=start, end=end, freq='1h')

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
P_PV_STC = 14.4                                                                 # rated power from TüV Calibration report under STC Conditions at 25°C & 1000 W/m^2                            
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

dwd_data = pd.read_excel(os.path.join("Imports", r'Timeseries_JRC_PVGIS_TH_Koeln.xlsx'))  # Hourly Weather Data (DNI , GHI , DHI , temp_air , wind speed and pressure)

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
house_data_read = pd.read_excel(os.path.join("Imports", r'house_demand.xlsx'))

house_data = pd.DataFrame(index=house_data_read.index, columns=["elec_cons", "thermal_cons"])

"_____________Assign Headers to pv_data Dataframe______________"
house_data["DateTimeIndex"] = house_data_read.date
house_data["DateTimeIndex"] = pd.to_datetime(house_data["DateTimeIndex"])
house_data["elec_cons"] = house_data_read.elec_cons
house_data["thermal_cons"] = house_data_read.thermal_cons

"______Import_Operating_strategies & Mass_Flow_Loss_values____________"

#Operating Strategies
op_strategy = pd.read_excel(os.path.join("Imports", r'rauhigkeit.xlsx'), sheet_name='Parameter')

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
    m_flow_loss_temp = pd.read_excel(os.path.join("Imports", r'rauhigkeit.xlsx'), sheet_name=os_name.format(str(i)), usecols = "A,I:J, L, R")
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
ks_SRT = 0.000225                                                               # ks/roughness value for one SRT, used in design mode to calculate the pressure drop. ks_SRT values for off design mode are calculated
p_amb=1.01325                                                                   # Atmospheric pressure [Bar]
mass_flow = 0.0320168                                                         # Can be one value or string (from measurement data later on). IMPORTANT: This mass flow value applies for one String of 12 SRTs and is not the mass flow delivered by the fan for the whole SRT plant!
P_HP = 700                                                                      #nominal power of the heat Pump in Watts, taken as constant

# Allowed Value range for mass flow is:
# 0.0646 to 0.00306 kg/s (due to interpolation boundaries)

mass_flow_temp = None                                                           #Can be ignored
E_sdp_Cooling = 0

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

"_____________Calculations____________"
"_____________Electrical Yield____________"

######
# Electrical Yield
######

dfMainElec = [] #overall electrical results
dfMainElecNew = [] # overall results with cooling effect
dfSubElec = [] # Row result
dfSubElec_New = [] # Row Result with cooling effect
dfThermalMain = [] #Thermal Results
dfThermalSub = [] # Thermal Effect of one row
df_test_elecindicators_full = [] # testing performance indicators
totalPowerDiff = 0

#for i in tqdm(pv_data.index[8:10]):
for i in tqdm(pv_data.index[5548:5564]):


    "_______Looping through excel rows_______"
    "Aligning excel row values to variable"

    time = pv_data.DateTimeIndex[i]
    temp_amb = pv_data.temp_air[i]
    wind_amb = pv_data.wind_speed[i]
    ghi = pv_data.ghi[i]
    dni = pv_data.dni[i]
    dhi = pv_data.dhi[i]
    Tamb = pv_data.temp_air[i]

    "______Calculating ks value & mass flow leakage via interpolation______"

    if mass_flow_temp != mass_flow:
        m_loss_offdesign, ks_SRT, mass_flow_temp = sdp.interpolate_ks_mloss(i=i,
                                                              op_strategy=op_strategy,
                                                              os_name=os_name,
                                                              mass_flow=mass_flow,
                                                              mass_flow_loss=mass_flow_loss,
                                                              mass_flow_temp=mass_flow_temp)
        m_loss_offdesign = m_loss_offdesign.set_axis(['5_1_dpx'], axis=1, inplace=False)


    "______Getting the initial cell temperature______"
    initCellTemperature = cellTemperature(latitude=latitude, longitude=longitude,
                                          m_azimut=m_azimut, m_tilt=m_tilt,
                                          time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
                                          dhi=pv_data.dhi[i],
                                          albedo=albedo, irrad_model=irrad_model,
                                          wind_amb=pv_data.wind_speed[i],
                                          temp_avg=pv_data.temp_air[i])


    "________Finding the electrical yield & Solar data based in initial cell temperature_______"
    electrical_yield = Photovoltaic(latitude=latitude, longitude=longitude, altitude=altitude, timezone=timezone,
                                    m_azimut=m_azimut, m_tilt=m_tilt, module_number=num_sdp_series*num_sdp_parallel,
                                    time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i], dhi=pv_data.dhi[i],
                                    albedo=albedo, a_r=a_r, irrad_model=irrad_model, module=module,
                                    temp_amb=pv_data.temp_air[i], wind_amb=pv_data.wind_speed[i],pressure=pv_data.pressure[i],cell_temp=initCellTemperature.tcell)

    "_______Making an Array of results got from electrical_yieldcmd_______"
    if int(electrical_yield.effective_irradiance) == 0 :
        efficency = 0
    else :
        efficency = round((round(electrical_yield.annual_energy, 2) / (
                num_sdp_series * num_sdp_parallel * 0.1)) * 100 / int(electrical_yield.effective_irradiance),2)
    
    "_____Calculating_Performance_Indicators_(electrical)_Values_with_Cooling_effect_________"
    # Calculation according to DIN EN IEC 61724-1 (VDE 0126-25-1): 2021-03 exept of the autarky rate
        
    if int(electrical_yield.effective_irradiance) == 0:
        E_ideal_STC = 0
        C_k_STC = 0
        PR_STC = 0
        C_k_annual = 0
        PR_annual_eq = 0
        Y_F = 0
        autarky_rate = 0
    else:
        # Bemessungsleitungs Temperaturanpassungsfaktor Ck
        # with y (module["gamma_pmp"]) as relativer Höchstleistungs-Temperaturkoeffizient in 1/C from TüV Report
        C_k_STC = 1 + (module["gamma_pmp"]/100) * (electrical_yield.tcell.item()-25)                                                    # for STC conditions
        C_k_annual = 1 + (module["gamma_pmp"]/100) * (electrical_yield.tcell.item()- 22)                                                # temperature corrected                                                                                                                                                                                   
        
        # Ideal energy yield under STC conditions and temperature corrected
        E_ideal_STC = round((((P_PV_STC * C_k_STC * (num_sdp_series * num_sdp_parallel)) * float(electrical_yield.effective_irradiance)) / 1000), 2)  
        
        # Perforance ratio under STC conditions and temperature corrected
        PR_STC = round((electrical_yield.power_ac / E_ideal_STC), 2)
        PR_annual_eq = round(electrical_yield.power_ac / (((P_PV_STC * C_k_annual * (num_sdp_series * num_sdp_parallel)) * float(electrical_yield.effective_irradiance)) / 1000), 2)
        
        # Final system yield (net energy yield)
        Y_F = round((float(electrical_yield.power_ac) / (P_PV_STC * (num_sdp_series * num_sdp_parallel))), 2)
                
        # Autarky rate according to Quaschning, V.: Regenerative Energiesysteme (2019)
        autarky_rate = round((float(electrical_yield.power_ac) / house_data["elec_cons"][i])*100, 2)
                    
        #df_test_elecindicators = [i, time, round(t_avg_new, 2), round(electrical_yield_new.annual_energy, 2), int(electrical_yield_new.effective_irradiance), E_ideal, PR, Y_F, autarkie_rate]
        
    # use t ambient for no cooling effect
    dfSubElec = [i, time, temp_amb, round(electrical_yield.annual_energy, 2),
                 int(electrical_yield.effective_irradiance),efficency,E_ideal_STC, PR_STC, PR_annual_eq, Y_F, autarky_rate,0,0]
    "_______Electrical Yield for one cell_____"
    P_MP = dfSubElec[3] / (num_sdp_series * num_sdp_parallel) #Anual Energy / total modules
    effective_Iradiance = dfSubElec[4]
    "______Overall per unit area____"
    E_sdp_New = (0.93 * (effective_Iradiance * module['Area']) - (P_MP)) / module['Area']

    if E_sdp_New == 0:
        # in deg Celsius
        t_out_init = Tamb

        # Watt
        p_fan_init = 0

        # kg/sec
        m_out_init = 0

    else:
        "_______SDP Calculations to find t_out, fan in , mass flow out"
        t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(
            ambient_temp=pv_data.temp_air[i],
            absorption_incl=E_sdp_New,
            inlet_temp=pv_data.temp_air[i],
            mass_flow=mass_flow,
            print_res=False,
            ks_SRT=ks_SRT,
            m_loss_offdesign=m_loss_offdesign,
            )

    t_out = t_out_init
    "____Finding Avg temperature for cooling effect_____"

    T_PV_Temp_Model = float(initCellTemperature.tcell)

    #  t_avg = (T_PV_Temp_Model + t_out) / 2

    t_m = (temp_amb + t_out) / 2

    t_avg_new = (T_PV_Temp_Model + t_m) / 2

    Residue =  T_PV_Temp_Model - t_avg_new

    totalLoops = 0

    while totalLoops < 1:
        totalLoops = totalLoops + 1
        newCellTemperature = cellTemperature(latitude=latitude, longitude=longitude,
                                             m_azimut=m_azimut, m_tilt=m_tilt,
                                             time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
                                             dhi=pv_data.dhi[i],
                                             albedo=albedo, irrad_model=irrad_model,
                                             wind_amb=pv_data.wind_speed[i],
                                             temp_avg=t_avg_new)

        t_avg_old = t_avg_new


        "________This electrical_yield_new is using cell temperature found from cooling effect above________"
        electrical_yield_new = Photovoltaic(latitude=latitude, longitude=longitude, altitude=altitude,
                                            timezone=timezone,
                                            m_azimut=m_azimut, m_tilt=m_tilt,
                                            module_number=num_sdp_series * num_sdp_parallel,
                                            time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
                                            dhi=pv_data.dhi[i],
                                            albedo=albedo, a_r=a_r, irrad_model=irrad_model, module=module,
                                            temp_amb=pv_data.temp_air[i], wind_amb=pv_data.wind_speed[i],
                                            pressure=pv_data.pressure[i], cell_temp=t_avg_new)
        "________Calculating the Heat Pump COP_with Cooling effect_______"
        heatPump = HeatPump(0.068)
        heatPumpCOP = round(heatPump.calc_cop_ruhnau(50, t_out_init, "ashp"), 2)
        heatPumpCOP_ref = round(heatPump.calc_cop_ruhnau(50, pv_data.temp_air[i], "ashp"), 2)

        "_______Calculating Heat Pump thermal Output_______"
        
        heatPumpThermal = round(68 * heatPumpCOP, 2)
        heatPumpThermal_ref = round(68 * heatPumpCOP_ref, 2)

        if int(electrical_yield_new.effective_irradiance) == 0:
            efficency_New = 0
        else:
            efficency_New = round((round(electrical_yield_new.annual_energy, 2)/(num_sdp_series * num_sdp_parallel*0.1))*100/int(electrical_yield_new.effective_irradiance),2)
            
        "_____Calculating_Performance_Indicators_(electrical)_Values_with_Cooling_effect_________"
        # Calculation according to DIN EN IEC 61724-1 (VDE 0126-25-1): 2021-03 exept of the autarky rate
        
        if int(electrical_yield_new.effective_irradiance) == 0:
            E_ideal_STC = 0
            C_k_STC = 0
            PR_STC = 0
            C_k_annual = 0
            PR_annual_eq = 0
            Y_F = 0
            autarky_rate = 0
        else:
            # Bemessungsleitungs Temperaturanpassungsfaktor Ck
            # with y (module["gamma_pmp"]) as relativer Höchstleistungs-Temperaturkoeffizient in 1/C from TüV Report
            C_k_STC = 1 + (module["gamma_pmp"]/100) * (t_avg_new-25)                                                    # for STC conditions
            C_k_annual = 1 + (module["gamma_pmp"]/100) * (t_avg_new- 22)                                                # temperature corrected                                                                                                                                                                                   
            
            # Ideal energy yield under STC conditions and temperature corrected
            E_ideal_STC = round((((P_PV_STC * C_k_STC * (num_sdp_series * num_sdp_parallel)) * float(electrical_yield_new.effective_irradiance)) / 1000), 2)  
            
            # Perforance ratio under STC conditions and temperature corrected
            PR_STC = round((electrical_yield_new.power_ac / E_ideal_STC), 2)
            PR_annual_eq = round(electrical_yield_new.power_ac / (((P_PV_STC * C_k_annual * (num_sdp_series * num_sdp_parallel)) * float(electrical_yield_new.effective_irradiance)) / 1000), 2)
            
            # Final system yield (net energy yield)
            Y_F = round((float(electrical_yield_new.power_ac) / (P_PV_STC * (num_sdp_series * num_sdp_parallel))), 2)
                    
            # Autarky rate according to Quaschning, V.: Regenerative Energiesysteme (2019)
            autarky_rate = round((float(electrical_yield_new.power_ac) / house_data["elec_cons"][i])*100, 2)
                        
        #df_test_elecindicators = [i, time, round(t_avg_new, 2), round(electrical_yield_new.annual_energy, 2), int(electrical_yield_new.effective_irradiance), C_k_STC, E_ideal_STC, PR_STC, C_k_annual, Y_F, PR_annual_eq, autarky_rate]
        
        "______Saving results__________"                          
        dfSubElec_New = [i, time, t_avg_new, round(electrical_yield_new.annual_energy, 2), int(electrical_yield_new.effective_irradiance),
                         efficency_New, E_ideal_STC, PR_STC, PR_annual_eq, Y_F, autarky_rate,heatPumpCOP,heatPumpThermal]

        P_MP_New = dfSubElec_New[3] / (num_sdp_series * num_sdp_parallel)
        effective_Iradiance_New = dfSubElec_New[4]

        # It will be lower for cooling
        E_sdp_Cooling = (0.93 * (effective_Iradiance_New * module['Area']) - (P_MP_New)) / module['Area']

        if E_sdp_Cooling == 0:
            # in deg Celsius
            t_out = Tamb

            # Watt
            p_fan = 0

            # kg/sec
            m_out = 0

        else:

            t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(
                ambient_temp=pv_data.temp_air[i],
                absorption_incl=E_sdp_Cooling,
                inlet_temp=pv_data.temp_air[i],
                mass_flow=mass_flow,
                print_res=False,
                ks_SRT=ks_SRT,
                m_loss_offdesign=m_loss_offdesign,
            )

        t_m = (temp_amb + t_out_init) / 2
        t_avg_new = (T_PV_Temp_Model + t_m) / 2

        Residue =  t_avg_old - t_avg_new


    "________Find power diffeence________"
    powerDiff = P_MP_New - P_MP
    totalPowerDiff = powerDiff + totalPowerDiff

    p_fan = p_fan_init
    m_out = m_out_init

    # Calculating the heat flux normed on one m^2 (Division through number of SRTs and their area 0.10)
    flux = round((m_out * 1.005 * (t_out - Tamb) / (num_sdp_series * num_sdp_parallel * 0.10)), 3) # cp_air: 1.005 kJ/kg*K, Unit is kJ/s --> kW
    
    "_____Calculating_Performance_Indicators_(thermal)_Values_only_with_Cooling_effect_________"
     
    #solarer_Deckungsgrad - without taking the storage into consideraiton
    #sol_coverage = round(((house_data["thermal_cons"][i] - heatPumpThermal)/house_data["thermal_cons"][i])*100, 2)
               
    elec_parameter = (house_data.elec_cons[i] + p_fan) < (P_MP_New * num_sdp_parallel * num_sdp_series)
    thermal_parameter = (house_data.thermal_cons[i] < flux)

    #E_sdp_Cooling = E_sdp_Cooling

    if E_sdp_Cooling == 0:
        status = "System Off"

    elif (elec_parameter is True) and (thermal_parameter is True):
        status = "Good"

    elif (elec_parameter is True) and (thermal_parameter is False):
        status = "Increase Thermal Production (increase fan_power)"


    elif (elec_parameter is False) and (thermal_parameter is True):
        status = "Decrease Thermal Production (reduce fan_power)"

    else:
        status = "System Off"

    if dfSubElec[4] == 0:
        thermal_efficency = 0
    else:
        thermal_efficency = round((flux * 1000*100) / dfSubElec_New[4],2)
    dfThermalSub = [i, time, Tamb, round(E_sdp_Cooling, 2), t_out, p_fan, m_out, flux,thermal_efficency, status, elec_parameter,
                    thermal_parameter]
    dfThermalMain.append(dfThermalSub)

    dfMainElec.append(dfSubElec)
    dfMainElecNew.append(dfSubElec_New)
    
    #df_test_elecindicators_full.append(df_test_elecindicators)
    
#df_test_elecindicators_full = pd.DataFrame(df_test_elecindicators_full)     

"_________Saving results in excel_____________"

column_values_elec = ["Index", "Time", "Tavg [°C]", "Power [W]", "Effective Irradiance [W/m^2]", "Elec. Efficency [%]","ideal elec. energy yield [Wh]","Performance Ratio STC[-]","Performance Ratio eq [-]","spez. elec. energy yield [kWh/kWp]","Autarky rate [%]","HP_COP [-]","HP_Thermal[Wh]"]
# Assigning df all data to new varaible electrical data
electrical_data_New = pd.DataFrame(data=dfMainElecNew, columns=column_values_elec)
electrical_data_New.fillna(0)  # fill empty rows with 0
div = len(electrical_data_New[electrical_data_New["Effective Irradiance [W/m^2]"] > 0])
electrical_data_New.loc['Total'] = electrical_data_New.select_dtypes(np.number).sum()  # finding total number of rows

"""
electrical_data_New.loc['Averages'] = electrical_data_New.select_dtypes(np.number).sum()
electrical_data_New.loc['Averages'] = electrical_data_New.divide(electrical_data_New, axis={"Effective Irradiance [W/m^2]", 
                     "Performance Ratio STC[-]", 
                     "Performance Ratio eq [-]",
                     "spez. elec. energy yield [kWh/kWp]",
                     "Autarky rate [%]",
                     "HP_COP [-]",
                     "HP_Thermal[Wh]"})

electrical_data_New['Averages'] = electrical_data_New.loc[['Averages'], ['Effective Irradiance [W/m^2]']].div(div)
electrical_data_New.loc['Averages'] = electrical_data_New.div(div, level=3)
"""

pd.set_option('display.max_colwidth', 40)
electrical_data_New.to_excel(os.path.join("Exports", r'ResultsWithCoolingEffect.xlsx'))


electrical_data = pd.DataFrame(data=dfMainElec, columns=column_values_elec)
electrical_data.fillna(0)  # fill empty rows with 0
electrical_data.loc['Total'] = electrical_data.select_dtypes(np.number).sum()  # finding total number of rows
#electrical_data.loc['Total']['Time'] = "Sum"
pd.set_option('display.max_colwidth', 40)
electrical_data.to_excel(os.path.join("Exports", r'ResultsWithoutCoolingEffect.xlsx'))


column_values = ["Index", "Time", "Tamb", "E_sdp_eff", "T_out", "P_fan", "M_out", "HeatFlux_[kW/m^2]","Thermal Efficency","status",
                 "Elec_demand_met", "Heat_demand_met"]
thermal_data = pd.DataFrame(data=dfThermalMain, columns=column_values)
thermal_data.loc['Total'] = thermal_data.select_dtypes(np.number).sum()
pd.set_option('display.max_colwidth', 8)
# print(thermal_data)

Efficiency = thermal_data.loc["Total", "HeatFlux_[kW/m^2]"] / thermal_data.loc["Total", "E_sdp_eff"]
print(f'Total Power difference with and without cooling effect {round(totalPowerDiff, 2)} Watt hours')
print("Efficiency wrt Effective Irradiance:", round(Efficiency * 100, 2), "%")
complete_data = pd.merge(electrical_data_New, thermal_data)
complete_data.to_excel(os.path.join("Exports", r'CompleteResult.xlsx'))

complete_data = pd.merge(electrical_data, thermal_data)
complete_data.to_excel(os.path.join("Exports", r'CompleteResultWithoutCoolingEffect.xlsx'))


"____Calculating_Jahresarbeitszahl_of the heat pump____"
Jahresarbeitszahl_ref = round(electrical_data_New.select_dtypes(np.number).sum().loc['HP_Thermal[Wh]'] / (P_HP * len(electrical_data_New)), 2)
Jahresarbeitszahl = round(electrical_data_New.select_dtypes(np.number).sum().loc['HP_Thermal[Wh]'] / (P_HP * len(electrical_data_New)), 2) 
print("\nAnnual performance factor (Jahresarbeitzahl) of the heat Pump: ", Jahresarbeitszahl, "Annual Performance factor reference of the heat Pump with Tamb: ", Jahresarbeitszahl_ref)


"________Plotting______"

P_Valve = sdp.plot_temperature_curve(p_amb=p_amb)

#len(electrical_data_New[electrical_data_New["Effective Irradiance [W/m^2]"] > 0])
