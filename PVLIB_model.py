# -*- coding: utf-8 -*-

import pvlib


# from pvlib.pvsystem import pvsystem
# pip install NREL-PySAM
# NREL PySAM has to be installed to run fit_cec_sam


class Photovoltaic():
    def __init__(self, latitude, longitude, altitude, timezone, m_azimut, m_tilt, module_number, time, dni, ghi, dhi,
                 albedo, a_r, irrad_model, module, temp_amb, wind_amb, pressure):
        # Get Parameters
        self.temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm'][
            'close_mount_glass_glass']

        self.times = time

        self.latitude = latitude
        self.longitude = longitude

        self.module_number = module_number
        self.m_azimut = m_azimut
        self.m_tilt = m_tilt

        self.clearsky_dni = dni
        self.clearsky_ghi = ghi
        self.clearsky_dhi = dhi

        self.albedo = albedo
        self.a_r = a_r
        self.irrad_model = irrad_model
        self.module = module

        # Radiation & Power Calculations
        self.solpos = pvlib.solarposition.get_solarposition(self.times, self.latitude, self.longitude)

        self.dni_extra = pvlib.irradiance.get_extra_radiation(self.times)

        self.airmass = pvlib.atmosphere.get_relative_airmass(self.solpos['apparent_zenith'])

        self.am_abs = pvlib.atmosphere.get_absolute_airmass(self.airmass, pressure)

        self.aoi = pvlib.irradiance.aoi(self.m_tilt, self.m_azimut,
                                        self.solpos['apparent_zenith'], self.solpos['azimuth'])

        self.total_irrad = pvlib.irradiance.get_total_irradiance(self.m_tilt, self.m_azimut,
                                                                 self.solpos['apparent_zenith'],
                                                                 self.solpos['azimuth'],
                                                                 self.clearsky_dni, self.clearsky_ghi,
                                                                 self.clearsky_dhi,
                                                                 dni_extra=self.dni_extra,
                                                                 model=self.irrad_model,
                                                                 albedo=self.albedo)

        # Estimating cell temperatue via Faimann Model
        self.tcell = pvlib.temperature.faiman(self.total_irrad['poa_global'], temp_amb, wind_amb, u0=25.0, u1=6.84)
        # Manual Effective Irradiances calculation via:
        "Ee = f_1(AM_a) (E_b f_2(AOI) + f_d E_d)"
        "or"
        "Ee = F1 * (poa_direct * F2 + module['FD'] * poa_diffuse)"

        self.IAM_mod_diff = pvlib.iam.martin_ruiz_diffuse(self.m_tilt, a_r=self.a_r, c1=0.4244,
                                                          c2=None)  # Estimating Angle of incidence modifier(IAM) using the Martin and Ruiz for diffuse radiation

        self.F1 = pvlib.pvsystem.sapm_spectral_loss(self.am_abs,
                                                    self.module)  # Calculating SAPM spectral loss coefficient, F1.

        self.IAM_mod_dir = pvlib.iam.martin_ruiz(self.aoi,
                                                 a_r=self.a_r)  # Estimating Angle of incidence modifier(IAM) using the Martin and Ruiz for direct radiation

        self.effective_irradiance = self.F1 * (self.total_irrad['poa_direct'] * self.IAM_mod_dir
                                               + self.IAM_mod_diff[0] * self.total_irrad['poa_sky_diffuse']
                                               + self.IAM_mod_diff[1] * self.total_irrad[
                                                   'poa_ground_diffuse'])  # Effective Irradiance calculation --> Ground reflection considered

        # Estimates parameters for the CEC single diode model (SDM) using the SAM SDK
        self.mp_fit_cec_sam = pvlib.ivtools.sdm.fit_cec_sam(
            self.module['celltype'],
            v_mp=self.module['Vmpo'],
            i_mp=self.module['Impo'],
            v_oc=self.module['Voco'],
            i_sc=self.module['Isco'],
            alpha_sc=self.module['Aisc'],
            beta_voc=self.module['Bvoco'],
            gamma_pmp=self.module['gamma_pmp'],
            cells_in_series=self.module['Cells_in_Series'],
            temp_ref=25)

        # Calculating parameters for solving the single diod model by De Soto et. al model
        self.single_diod_parameters = pvlib.pvsystem.calcparams_desoto(
            self.effective_irradiance,
            self.tcell,
            alpha_sc=self.module['Aisc'],
            a_ref=self.mp_fit_cec_sam[4],
            I_L_ref=self.mp_fit_cec_sam[0],
            I_o_ref=self.mp_fit_cec_sam[1],
            R_sh_ref=self.mp_fit_cec_sam[3],
            R_s=self.mp_fit_cec_sam[2],
            EgRef=1.121, dEgdT=- 0.0002677, irrad_ref=1000, temp_ref=25)

        # Solving the single-diode equation to obtain a photovoltaic IV curve. Important is the output p_mp
        self.single_diod_IV = pvlib.pvsystem.singlediode(
            photocurrent=self.single_diod_parameters[0],
            saturation_current=self.single_diod_parameters[1],
            resistance_series=self.single_diod_parameters[2],
            resistance_shunt=self.single_diod_parameters[3],
            nNsVth=self.single_diod_parameters[4],
            ivcurve_pnts=None,
            method='lambertw')

        # Storing module Power to pass it on to the Excel Sheet
        self.annual_energy = self.module_number * self.single_diod_IV[
            'p_mp'].sum()  # in Watt hours, because p_mp is given in Watt



class cellTemperature():
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

