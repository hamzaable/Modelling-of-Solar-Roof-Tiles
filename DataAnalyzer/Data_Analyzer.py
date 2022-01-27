import os

import pandas as pd
import numpy as np

from tqdm import tqdm

#
filename = "ToAnalyze.xlsx"
column_values = ["Index", "Time", "Tamb [°C]","Tmcooling [°C]","Tm [°C]","T heatflux [°C]", "Power-DC [W]", "Power-AC [W]", "Effective Irradiance [W/m^2]", "Elec. Efficency [%]","ideal elec. energy yield [Wh]","Performance Ratio [-]","Performance Ratio STC[-]","Performance Ratio eq [-]","spez. elec. energy yield [kWh/kWp]","Autarky rate [%]","HP_COP [-]","HP_Thermal[Wh]"]
resultSheetName = "NewResults.xlsx"

#

plain_data = pd.read_excel(filename)
data = pd.DataFrame(index=plain_data.index,data=plain_data, columns=column_values)

for i in tqdm(data.index[0:100]): #8759
    ambiant_temperature = data["Tamb [°C]"][i]*100
    data["Tamb [°C]"][i] = ambiant_temperature
    heat_flux = data.Time[i]


data.to_excel(resultSheetName)