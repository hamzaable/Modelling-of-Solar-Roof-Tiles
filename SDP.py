# -*- coding: utf-8 -*-
#  Test cooling effect
import pandas as pd
import numpy as np

##########
from PVLIB_model import Photovoltaic
from TESPy_model import SDP_sucking

##########
start = "01-01-{} 00:00".format(str(2019))
end = "31-12-{} 23:00".format(str(2019))
naive_times = pd.date_range(start=start, end=end, freq='1h')

latitude = 50.9375
longitude = 6.9603
name = 'Cologne'
altitude = 121
timezone = 'Etc/GMT+2'

dwd_data = pd.read_excel(r'704EEE00.xlsx')  # Hourly Weather Data (DNI , GHI , DHI , temp_air , wind speed and pressure)

pv_data = pd.DataFrame(index=dwd_data.index, columns=["dni", "ghi",
                                                      "dhi",
                                                      "temp_air",
                                                      "wind_speed", "pressure"])
# Assign Headers to pv_data Dataframe

pv_data["DateTimeIndex"] = dwd_data.date
pv_data["DateTimeIndex"] = pd.to_datetime(pv_data["DateTimeIndex"])
pv_data["dni"] = dwd_data.irradiance_dir
pv_data["dhi"] = dwd_data.irradiance_diff
pv_data["ghi"] = dwd_data.poa_global
pv_data["temp_air"] = dwd_data.temp
pv_data["wind_speed"] = dwd_data.wind_speed
pv_data["pressure"] = dwd_data.pr


house_data_read = pd.read_excel(r'house_demand.xlsx') # Hourly house demand (elec_cons , thermal_cons)
house_data = pd.DataFrame(index=house_data_read.index, columns=["elec_cons", "thermal_cons"])

# Assign Headers to pv_data Dataframe
house_data["DateTimeIndex"] = house_data_read.date
house_data["DateTimeIndex"] = pd.to_datetime(house_data["DateTimeIndex"])
house_data["elec_cons"] = house_data_read.elec_cons
house_data["thermal_cons"] = house_data_read.thermal_cons

num_sdp_series = 12
num_sdp_parallel = 38

#####
# Thermal initialization
#####
sdp = SDP_sucking(sdp_in_parallel=num_sdp_parallel,
                  sdp_in_series=num_sdp_series)

sdp.init_sdp(ambient_temp=-4,
             absorption_incl=300,
             inlet_temp=-4,
             mass_flow=1,
             # zeta=4e6,
             m_loss=0.001,
             print_res=False)

######
# Electrical Yeild
######

df = []
df1 = []

dfMainElec = []
dfSubElec = []
dfThermalMain = []
dfThermalSub = []
# for i in pv_data.index[0:8760]:
# Looping through weather data profile
for i in pv_data.index[0:30]:
    time = pv_data.DateTimeIndex[i]
    temp_amb = pv_data.temp_air[i]
    wind_amb = pv_data.wind_speed[i]
    ghi = pv_data.ghi[i]
    dni = pv_data.dni[i]
    dhi = pv_data.dhi[i]
    Tamb = pv_data.temp_air[i]
    # finding different weather parameters using pvlib
    electrical_yield = Photovoltaic(latitude=latitude, longitude=longitude, altitude=altitude, timezone=timezone,
                                    time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
                                    dhi=pv_data.dhi[i], temp_amb=pv_data.temp_air[i], wind_amb=pv_data.wind_speed[i],
                                    pressure=pv_data.pressure[i])
    # Making an Array of results got from electrical_yield
    dfSubElec = [i, time, temp_amb, round(electrical_yield.annual_energy, 2), int(electrical_yield.effective_irradiance)]

    # This point is important because here we can add cooling effect
    E_sdp_1 = (0.93 * (dfSubElec[4] * 0.10)
               - (dfSubElec[3])) / 0.10

    if E_sdp_1 == 0:
        # in deg celcius
        t_out_init = Tamb

        # Watt
        p_fan_init = 0

        # kg/sec
        m_out_init = 0

    else:
        t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(
            ambient_temp=pv_data.temp_air[i],
            absorption_incl=E_sdp_1,
            inlet_temp=pv_data.temp_air[i],
            mass_flow=1,
            print_res=False)

    t_out = t_out_init
    p_fan = p_fan_init
    m_out = m_out_init
    flux = round((m_out * (t_out - Tamb) / (num_sdp_series * num_sdp_parallel * 0.10)), 2)

    elec_parameter = (house_data.elec_cons[i] + p_fan) \
                     < (dfSubElec[3] * num_sdp_parallel * num_sdp_series)
    thermal_parameter = (house_data.thermal_cons[i] < flux)

    if E_sdp_1 == 0:
        status = "System Off"

    elif (elec_parameter is True) and (thermal_parameter is True):
        status = "Good"

    elif (elec_parameter is True) and (thermal_parameter is False):
        status = "Increase Thermal Production (increase fan_power)"
        # t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(ambient_temp=pv_data.temp_air[i],
        # absorption_incl=pv_data.ghi[i], inlet_temp=pv_data.temp_air[i], mass_flow=0.42, print_res=False) t_out =
        # t_out_init p_fan =  p_fan_init m_out =  m_out_init flux = round((m_out*(t_out-Tamb)/(
        # num_sdp_series*num_sdp_parallel*0.10)),2)

        # elec_parameter = (house_data.elec_cons[i]+p_fan)<(electrical_data.Power[i]*num_sdp_parallel*num_sdp_series)
        # thermal_parameter = (house_data.thermal_cons[i]<flux)

        # elec_parameter = (house_data.elec_cons[i]+p_fan)<(electrical_data.Power[i]*num_sdp_parallel*num_sdp_series)
        # thermal_parameter = (house_data.thermal_cons[i]<flux)

    elif (elec_parameter is False) and (thermal_parameter is True):
        status = "Decrease Thermal Production (reduce fan_power)"
        # t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(ambient_temp=pv_data.temp_air[i],
        # absorption_incl=pv_data.ghi[i], inlet_temp=pv_data.temp_air[i], mass_flow=0.42, print_res=False)
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

    dfThermalSub = [i, time, Tamb, round(E_sdp_1, 2), t_out, p_fan, m_out, flux, status, elec_parameter, thermal_parameter]
    dfThermalMain.append(dfThermalSub)

    dfMainElec.append(dfSubElec)

column_values_elec = ["Index", "Time", "Tamb", "Power", "Effective_Irradiance"]
# Assigning df all data to new varaible electrical data
electrical_data = pd.DataFrame(data=dfMainElec, columns=column_values_elec)
electrical_data.fillna(0) # fill empty rows with 0
electrical_data.loc['Total'] = electrical_data.select_dtypes(np.number).sum() #  finding total number of rows
pd.set_option('display.max_colwidth', 40)
print(electrical_data)
electrical_data.to_excel(r'Results1.xlsx')

column_values = ["Index", "Time", "Tamb", "E_sdp_eff", "T_out", "P_fan", "M_out", "HeatFlux", "status",
                 "Elec_demand_met", "Heat_demand_met"]
thermal_data = pd.DataFrame(data=dfThermalMain, columns=column_values)
thermal_data.loc['Total'] = thermal_data.select_dtypes(np.number).sum()
pd.set_option('display.max_colwidth', 8)
print(thermal_data)

Efficiency = thermal_data.loc["Total", "HeatFlux"] / thermal_data.loc["Total", "E_sdp_eff"]
print("Efficiency wrt Effective Irradiance:", round(Efficiency * 100, 2), "%")

complete_data = pd.merge(electrical_data, thermal_data)
complete_data.to_excel(r'Result2.xlsx')

