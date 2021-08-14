import pandas as pd
import numpy as np


from tqdm import tqdm

##########
from PVLIB_model import Photovoltaic
from PVLIB_model import cellTemperature
from TESPy_model import SDP_sucking


class cooling_effect():

    # fint find the cell temperature
    def __init__(self, latitude, longitude, m_azimut, m_tilt, time, dni, ghi, dhi,
                 albedo, irrad_model, wind_amb, temp_avg):
        # Get Parameters


        self.times = time

        self.latitude = latitude
        self.longitude = longitude

        self.m_azimut = m_azimut
        self.m_tilt = m_tilt

        self.clearsky_dni = dni
        self.clearsky_ghi = ghi
        self.clearsky_dhi = dhi

        self.albedo = albedo
        self.irrad_model = irrad_model

        # Radiation & Power Calculations
        self.solpos = pvlib.solarposition.get_solarposition(self.times, self.latitude, self.longitude)

        self.dni_extra = pvlib.irradiance.get_extra_radiation(self.times)

        self.total_irrad = pvlib.irradiance.get_total_irradiance(self.m_tilt, self.m_azimut,
                                                                 self.solpos['apparent_zenith'],
                                                                 self.solpos['azimuth'],
                                                                 self.clearsky_dni, self.clearsky_ghi,
                                                                 self.clearsky_dhi,
                                                                 dni_extra=self.dni_extra,
                                                                 model=self.irrad_model,
                                                                 albedo=self.albedo)

        # Estimating cell temperatue via Faimann Model
        self.tcell = pvlib.temperature.faiman(self.total_irrad['poa_global'],
                                                  temp_avg,
                                                  wind_amb,
                                                  u0=25.0,
                                                  u1=6.84)



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

