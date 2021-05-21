# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 22:06:11 2020

@author: Vaishnavi Phadke
"""
import pvlib
#from pvlib.pvsystem import pvsystem
#pip install NREL-PySAM
# NREL PySAM has to be installed to run fit_cec_sam

#typical
# module = {"Vintage": 2020, "Area": 0.16, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.0010, "Aimp": -0.0003, "Bvoco": -0.0158, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "gamma_pmp": -0.3792, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}
# ar=0.25
# module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.0010, "Aimp": -0.0003, "Bvoco": -0.0158, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "gamma_pmp": -0.3792, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.0041, "B1": -0.002, 
#           "B2": 0.0002, "B3": -0.000008, "B4": 0.0000001, "B5":  -0.0000000009, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}
# #ar=0.18
# module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.0010, "Aimp": -0.0003, "Bvoco": -0.01608, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "gamma_pmp": -0.3792, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.0107, "B1": -0.00002, 
#           "B2":  0.0005, "B3": -2E-05, "B4": 3E-07, "B5": -2E-09, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}
# # =============================================================================
# #ar=0.14
module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8, 
          "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
          "Aisc": 0.0010, "Aimp": -0.0003, "Bvoco": -0.0158, "Mbvoc": 0, "Bvmpo": -0.01608, 
          "Mbvmp": 0, "gamma_pmp": -0.3792, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
          "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.0189, "B1": -0.0089, 
          "B2":0.0009 , "B3": -3E-05, "B4": 5E-07, "B5": -3E-09, "DTC": 3, "FD": 1, 
          "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
          "C4": 1, "C5": 0, "C6": 1, "C7": 0}
# =============================================================================
# #Alberino P
# module = {"Vintage": 2020, "Material": "mc-Si", "celltype": "monoSi", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.0010, "Aimp": -0.0003, "Bvoco": -0.0158, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "gamma_pmp": -0.3792, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.069, "B1": -0.0251, 
#           "B2": 0.0022, "B3": -8E-05, "B4": 1E-06, "B5": -6E-09, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}

#module_params = list(module.values())



class Photovoltaic():
         def __init__(self, latitude, longitude, altitude, timezone, m_azimut, m_tilt, module_number, time, dni, ghi, dhi, temp_amb, wind_amb, pressure):
            
            #Get Parameters
            self.temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['close_mount_glass_glass']
            
            self.times = time
            
            self.latitude = latitude
            self.longitude = longitude
            
            self.module_number = module_number
            self.m_azimut =  m_azimut
            self.m_tilt = m_tilt
            
            self.clearsky_dni = dni
            self.clearsky_ghi = ghi
            self.clearsky_dhi = dhi
                
            #Radiation & Power Calculations
            self.solpos = pvlib.solarposition.get_solarposition(self.times, self.latitude, self.longitude)
            
            self.dni_extra = pvlib.irradiance.get_extra_radiation(self.times)        
            
            self.airmass = pvlib.atmosphere.get_relative_airmass(self.solpos['apparent_zenith'])
                    
            self.am_abs = pvlib.atmosphere.get_absolute_airmass(self.airmass, pressure)
            
            self.aoi = pvlib.irradiance.aoi(self.m_tilt, self.m_azimut,
                                   self.solpos['apparent_zenith'], self.solpos['azimuth'])
    
            self.total_irrad = pvlib.irradiance.get_total_irradiance(self.m_tilt, self.m_azimut,
                                                            self.solpos['apparent_zenith'],
                                                            self.solpos['azimuth'],
                                                            self.clearsky_dni, self.clearsky_ghi, self.clearsky_dhi,
                                                            dni_extra=self.dni_extra,
                                                            model='haydavies')
        
                        
            #Estimating cell temperatue via Faimann Model
            self.tcell = pvlib.temperature.faiman(self.total_irrad['poa_global'], temp_amb, wind_amb, u0=25.0, u1=6.84)
            
            "self.effective_irradiance = pvlib.pvsystem.sapm_effective_irradiance("
            "self.total_irrad['poa_direct'], self.total_irrad['poa_diffuse'],"
            "self.am_abs, self.aoi, module)"
            
            
            #Manual Effective Irradiances calculation via: 
            "Ee = f_1(AM_a) (E_b f_2(AOI) + f_d E_d)"
            "or"
            "Ee = F1 * (poa_direct * F2 + module['FD'] * poa_diffuse)"
            
            self.IAM_mod_diff = pvlib.iam.martin_ruiz_diffuse(self.m_tilt, a_r=0.14, c1=0.4244, c2=None)        #Estimating Angle of incidence modifier(IAM) using the Martin and Ruiz for diffuse radiation
            
            self.F1 = pvlib.pvsystem.sapm_spectral_loss(self.am_abs, module)                                    #Calculating SAPM spectral loss coefficient, F1.
                       
            self.IAM_mod_dir = pvlib.iam.martin_ruiz(self.aoi, a_r=0.14)                                        #Estimating Angle of incidence modifier(IAM) using the Martin and Ruiz for direct radiation
            
            self.effective_irradiance  = self.F1 * (self.total_irrad['poa_direct'] * self.IAM_mod_dir 
                                                    + self.IAM_mod_diff[0] * self.total_irrad['poa_diffuse'])   #Effective Irradiance calculation
                        
                      
            
            #Estimates parameters for the CEC single diode model (SDM) using the SAM SDK
            self.mp_fit_cec_sam = pvlib.ivtools.sdm.fit_cec_sam(
                module['celltype'], 
                v_mp=module['Vmpo'], 
                i_mp=module['Impo'], 
                v_oc=module['Voco'], 
                i_sc=module['Isco'],
                alpha_sc=module['Aisc'], 
                beta_voc=module['Bvoco'], 
                gamma_pmp=module['gamma_pmp'], 
                cells_in_series=module['Cells_in_Series'],
                temp_ref=25)
            
                               
            #Calculating parameters for solving the single diod model by De Soto et. al model
            self.single_diod_parameters = pvlib.pvsystem.calcparams_desoto(
                self.effective_irradiance, 
                self.tcell, 
                alpha_sc=module['Aisc'], 
                a_ref=self.mp_fit_cec_sam[4], 
                I_L_ref=self.mp_fit_cec_sam[0],
                I_o_ref=self.mp_fit_cec_sam[1], 
                R_sh_ref=self.mp_fit_cec_sam[3],
                R_s=self.mp_fit_cec_sam[2], 
                EgRef=1.121, dEgdT=- 0.0002677, irrad_ref=1000, temp_ref=25)
                        
            
            #Solving the single-diode equation to obtain a photovoltaic IV curve. Important is the output p_mp
            self.single_diod_IV = pvlib.pvsystem.singlediode(
                photocurrent=self.single_diod_parameters[0], 
                saturation_current=self.single_diod_parameters[1],
                resistance_series=self.single_diod_parameters[2], 
                resistance_shunt=self.single_diod_parameters[3],
                nNsVth=self.single_diod_parameters[4], 
                ivcurve_pnts=None, 
                method='lambertw')
            
            #Storing module Power to pass it on to the Excel Sheet 
            self.annual_energy = self.module_number*self.single_diod_IV['p_mp'].sum() #in Watt hours, because p_mp is given in Watt
            
        
