import os

import pandas as pd
import numpy as np

from tqdm import tqdm

#
filename = "ToAnalyze.xlsx"
column_values = ["Index", "Time", "Tamb [°C]","Tmcooling [°C]","Tm [°C]","T heatflux [°C]", "Power-DC [W]", "Power-AC [W]", "Effective Irradiance [W/m^2]", "Elec. Efficency [%]","ideal elec. energy yield [Wh]","Performance Ratio [-]","Performance Ratio STC[-]","Performance Ratio eq [-]","spez. elec. energy yield [kWh/kWp]","Autarky rate [%]","HP_COP [-]","HP_Thermal[Wh]"]

resultSheetName = "NewResults.xlsx"

module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8,
          "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, "Aisc": 0.0010, "Bvoco": -0.0158,
          "gamma_pmp": -0.3792, "A0": 0.9645, "A1": 0.02753, "A2": -0.002848, "A3": -0.0001439,
          "A4": 0.00002219}

P_PV_STC = 14.4

num_sdp_series = 96  # Changed from 12 to 2 for test purpose
num_sdp_parallel = 2  # Changed from 38 to 1 for test purpose

#

plain_data = pd.read_excel(filename)
data = pd.DataFrame(index=plain_data.index,data=plain_data, columns=column_values)
totalTemperature = data["Tmcooling [°C]"][8759]
totalTemperatureWithoutCooling = data["Tm [°C]"][8759]

# Make New columns
data["Performance Ratio eq [-] Cooling 8760"] =  pd.NaT
data["Performance Ratio eq [-] NotCooling 8760"] =  pd.NaT

for i in tqdm(data.index[0:8759]): #8759

    t_cooling = data["Tmcooling [°C]"][i]
    power_ac = data["Power-AC [W]"][i]
    effective_irradiance = data["Effective Irradiance [W/m^2]"][i]

    C_k_annual_cooling = 1 + (module["gamma_pmp"] / 100) * (
            t_cooling - totalTemperature/8760)
    C_k_annual_withoutcooling = 1 + (module["gamma_pmp"] / 100) * (
            t_cooling - totalTemperatureWithoutCooling/8760)

    if float(effective_irradiance) == 0 :
        PR_annual_eq_cooling = 0
        PR_annual_eq_withoutcooling = 0
    else :
        PR_annual_eq_cooling = round(power_ac / (((P_PV_STC * C_k_annual_cooling * (
                num_sdp_series * num_sdp_parallel)) * float(effective_irradiance)) / 1000),
                             5)
        PR_annual_eq_withoutcooling = round(power_ac / (((P_PV_STC * C_k_annual_withoutcooling * (
                num_sdp_series * num_sdp_parallel)) * float(effective_irradiance)) / 1000),
                             5)

    data["Performance Ratio eq [-] Cooling 8760"][i] = PR_annual_eq_cooling
    data["Performance Ratio eq [-] NotCooling 8760"][i] = PR_annual_eq_withoutcooling



data.to_excel(resultSheetName)

print("File Exported")