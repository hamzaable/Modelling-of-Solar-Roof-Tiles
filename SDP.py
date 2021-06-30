# -*- coding: utf-8 -*-
#  Test cooling effect
import pandas as pd
import numpy as np


from tqdm import tqdm

##########
from PVLIB_model import Photovoltaic
from PVLIB_model import cellTemperature
from TESPy_model import SDP_sucking

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

albedo = 0.20  # Ground reflection albedo factor 0.20 (Beton) --> K. Mertens- Photovoltaik S.52
a_r = 0.14  # Spectral Corrections factor for different module glasses
irrad_model = 'haydavies'  # Model for Irradiation calculation. Choose from: 'isotropic', 'klucher', 'haydavies', 'reindl', 'king', 'perez'
m_azimut = 180  # Module Azimut (Ausrichtung) [°]dwd_data = pd.read_excel(r'704EEE00.xlsx')  # Hourly Weather Data (DNI , GHI , DHI , temp_air , wind speed and pressure)
m_tilt = 45  # Module tilt (Neigung) [°]

# ======Module Parameters=======================================================
# #ar=0.14
module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8,
          "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568,
          "Aisc": 0.0010, "Bvoco": -0.0158, "Bvmpo": -0.01608, "gamma_pmp": -0.3792,
          "A0": 0.9645, "A1": 0.02753, "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219}
# =============================================================================

"_____________Data Imports_____________"

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

house_data_read = pd.read_excel(r'house_demand.xlsx')  # Hourly house demand (elec_cons , thermal_cons)
house_data = pd.DataFrame(index=house_data_read.index, columns=["elec_cons", "thermal_cons"])

# Assign Headers to pv_data Dataframe
house_data["DateTimeIndex"] = house_data_read.date
house_data["DateTimeIndex"] = pd.to_datetime(house_data["DateTimeIndex"])
house_data["elec_cons"] = house_data_read.elec_cons
house_data["thermal_cons"] = house_data_read.thermal_cons

"______TESPy Model Parameters_________"
num_sdp_series = 12  # Changed from 12 to 2 for test purpose
num_sdp_parallel = 16  ##Changed from 38 to 1 for test purpose   => 16

#####
# Thermal initialization
#####
sdp = SDP_sucking(sdp_in_parallel=num_sdp_parallel,
                  sdp_in_series=num_sdp_series)

sdp.init_sdp(ambient_temp=-4,
             absorption_incl=300,
             inlet_temp=-4,
             mass_flow=1,
             # zeta=890,
             m_loss=0.001,
             print_res=False)

"_____________Calculations____________"

######
# Electrical Yeild
######

dfMainElec = []
dfMainElecNew = []
dfSubElec = []
dfSubElec_New = []
dfThermalMain = []
dfThermalSub = []
totalPowerDiff = 0



for i in tqdm(pv_data.index[8:60]):
# for i in tqdm(pv_data.index[1:8760]):
    time = pv_data.DateTimeIndex[i]
    temp_amb = pv_data.temp_air[i]
    wind_amb = pv_data.wind_speed[i]
    ghi = pv_data.ghi[i]
    dni = pv_data.dni[i]
    dhi = pv_data.dhi[i]
    Tamb = pv_data.temp_air[i]

# Cooling Effect Calculations Start from here
    initCellTemperature = cellTemperature(latitude=latitude, longitude=longitude,
                                          m_azimut=m_azimut, m_tilt=m_tilt,
                                          time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
                                          dhi=pv_data.dhi[i],
                                          albedo=albedo, irrad_model=irrad_model,
                                          wind_amb=pv_data.wind_speed[i],
                                          temp_avg=pv_data.temp_air[i])


    # This electrical_yield is just for comparing power with and without cooling effect otherwise it has no use
    electrical_yield = Photovoltaic(latitude=latitude, longitude=longitude, altitude=altitude, timezone=timezone,
                                    m_azimut=m_azimut, m_tilt=m_tilt, module_number=num_sdp_series * num_sdp_parallel,
                                    time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
                                    dhi=pv_data.dhi[i],
                                    albedo=albedo, a_r=a_r, irrad_model=irrad_model, module=module,
                                    temp_amb=pv_data.temp_air[i], wind_amb=pv_data.wind_speed[i],
                                    pressure=pv_data.pressure[i],cell_temp=initCellTemperature.tcell)

    # Making an Array of results got from electrical_yield
    dfSubElec = [i, time, temp_amb, round(electrical_yield.annual_energy, 2),
                 int(electrical_yield.effective_irradiance)]

    P_MP = dfSubElec[3] / (num_sdp_series * num_sdp_parallel)
    effective_Iradiance = dfSubElec[4]
    E_sdp_New = (0.93 * (effective_Iradiance * module['Area']) - (P_MP)) / module['Area']

    if E_sdp_New == 0:
        # in deg Celsius
        t_out_init = Tamb

        # Watt
        p_fan_init = 0

        # kg/sec
        m_out_init = 0

    else:
        t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(
            ambient_temp=pv_data.temp_air[i],
            absorption_incl=E_sdp_New,
            inlet_temp=pv_data.temp_air[i],
            mass_flow=1,
            print_res=False)



    # Step 3
    t_out = t_out_init
    # Step 4
    T_PV_Temp_Model = float(initCellTemperature.tcell)
    t_avg = (T_PV_Temp_Model + t_out) / 2
    # Step 5 residue to perform ite
    t_avg_new = t_avg
    Residue = t_avg_new - t_out


    # while Residue > 0.25:
    #     # totalLoops = totalLoops + 1
    #     newCellTemperature = cellTemperature(latitude=latitude, longitude=longitude,
    #                                          m_azimut=m_azimut, m_tilt=m_tilt,
    #                                          time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
    #                                          dhi=pv_data.dhi[i],
    #                                          albedo=albedo, irrad_model=irrad_model,
    #                                          wind_amb=pv_data.wind_speed[i],
    #                                          temp_avg=t_avg_new)
    #
    #     t_avg_old = t_avg_new
    #     t_avg_new = (float(newCellTemperature.tcell) + t_out) / 2
    #     Residue = t_avg_new - t_avg_old
    #     # print(f'new residual {round(Residue, 4)}')
    #     if Residue < 0.25:
    #         break

    # print("total small loops = {}".format(totalLoops))

    # This electrical_yield_new is using ambiant temp found from cooling effect above
    electrical_yield_new = Photovoltaic(latitude=latitude, longitude=longitude, altitude=altitude, timezone=timezone,
                                        m_azimut=m_azimut, m_tilt=m_tilt,
                                        module_number=num_sdp_series * num_sdp_parallel,
                                        time=pv_data.DateTimeIndex[i], dni=pv_data.dni[i], ghi=pv_data.ghi[i],
                                        dhi=pv_data.dhi[i],
                                        albedo=albedo, a_r=a_r, irrad_model=irrad_model, module=module,
                                        temp_amb=pv_data.temp_air[i], wind_amb=pv_data.wind_speed[i],
                                        pressure=pv_data.pressure[i],cell_temp=t_avg_new)

    dfSubElec_New = [i, time, t_avg_new, round(electrical_yield_new.annual_energy, 2),
                     int(electrical_yield_new.effective_irradiance)]

    P_MP_New = dfSubElec_New[3] / (num_sdp_series * num_sdp_parallel)
    effective_Iradiance_New = dfSubElec_New[4]
    E_sdp_Cooling = (0.93 * (effective_Iradiance_New * module['Area']) - (P_MP_New)) / module['Area']

    if E_sdp_Cooling == 0:
        # in deg Celsius
        t_out_init = Tamb

        # Watt
        p_fan_init = 0

        # kg/sec
        m_out_init = 0

    else:
        t_out_init, p_fan_init, m_out_init = sdp.calculate_sdp(
            ambient_temp=pv_data.temp_air[i],
            absorption_incl=E_sdp_Cooling,
            inlet_temp=pv_data.temp_air[i],
            mass_flow=1,
            print_res=False)



    powerDiff = P_MP_New - P_MP
    totalPowerDiff = powerDiff + totalPowerDiff
    # Step 6 New E_SDP
    # E_sdp_New = (0.93 * (effective_Iradiance * module['Area']) - (P_MP_New)) / module['Area']
    # End of all steps

    p_fan = p_fan_init
    m_out = m_out_init

    flux = round((m_out * (t_out - Tamb) / (num_sdp_series * num_sdp_parallel * 0.10)), 2)

    elec_parameter = (house_data.elec_cons[i] + p_fan) \
                     < (P_MP_New * num_sdp_parallel * num_sdp_series)
    thermal_parameter = (house_data.thermal_cons[i] < flux)

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

    dfThermalSub = [i, time, Tamb, round(E_sdp_Cooling, 2), t_out, p_fan, m_out, flux, status, elec_parameter,
                    thermal_parameter]
    dfThermalMain.append(dfThermalSub)

    dfMainElec.append(dfSubElec)
    dfMainElecNew.append(dfSubElec_New)

column_values_elec = ["Index", "Time", "Tamb [°C]", "Power [W]", "Effective Irradiance [W/m^2]"]
# Assigning df all data to new varaible electrical data
electrical_data_New = pd.DataFrame(data=dfMainElecNew, columns=column_values_elec)
electrical_data_New.fillna(0)  # fill empty rows with 0
electrical_data_New.loc['Total'] = electrical_data_New.select_dtypes(np.number).sum()  # finding total number of rows
pd.set_option('display.max_colwidth', 40)
electrical_data_New.to_excel(r'ResultsWithCoolingEffect.xlsx')

# For Without Cooling Effect
electrical_data = pd.DataFrame(data=dfMainElec, columns=column_values_elec)
electrical_data.fillna(0)  # fill empty rows with 0
electrical_data.loc['Total'] = electrical_data.select_dtypes(np.number).sum()  # finding total number of rows
pd.set_option('display.max_colwidth', 40)
electrical_data.to_excel(r'ResultsWithoutCoolingEffect.xlsx')


# electrical_data.to_excel(r'Electrical_Yield_Single_Diod_Model.xlsx')

column_values = ["Index", "Time", "Tamb", "E_sdp_eff", "T_out", "P_fan", "M_out", "HeatFlux", "status",
                 "Elec_demand_met", "Heat_demand_met"]
thermal_data = pd.DataFrame(data=dfThermalMain, columns=column_values)
thermal_data.loc['Total'] = thermal_data.select_dtypes(np.number).sum()
pd.set_option('display.max_colwidth', 8)
# print(thermal_data)

Efficiency = thermal_data.loc["Total", "HeatFlux"] / thermal_data.loc["Total", "E_sdp_eff"]
print(f'Total Power difference with and without cooling effect {round(totalPowerDiff, 2)} Watt hours')
print("Efficiency wrt Effective Irradiance:", round(Efficiency * 100, 2), "%")
complete_data = pd.merge(electrical_data_New, thermal_data)
complete_data.to_excel(r'CompleteResult.xlsx')
complete_data = pd.merge(electrical_data, thermal_data)
complete_data.to_excel(r'CompleteResultWithoutCoolingEffect.xlsx')