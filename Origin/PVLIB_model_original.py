# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 22:06:11 2020

@author: Vaishnavi Phadke
"""
import pvlib
#from pvlib.pvsystem import pvsystem


latitude = 50.934055
longitude = 6.990349
name='Cologne'
altitude=  121
timezone= 'Etc/GMT+2'

#typical
# module = {"Vintage": 2020, "Area": 0.16, "Material": "mc-Si", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.00279, "Aimp": -0.0003, "Bvoco": -0.01608, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}
# ar=0.25
# module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.00279, "Aimp": -0.0003, "Bvoco": -0.01608, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.0041, "B1": -0.002, 
#           "B2": 0.0002, "B3": -0.000008, "B4": 0.0000001, "B5":  -0.0000000009, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}
# #ar=0.18
# module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.00279, "Aimp": -0.0003, "Bvoco": -0.01608, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.0107, "B1": -0.00002, 
#           "B2":  0.0005, "B3": -2E-05, "B4": 3E-07, "B5": -2E-09, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}
# # =============================================================================
# #ar=0.14
module = {"Vintage": 2020, "Area": 0.1, "Material": "mc-Si", "Cells_in_Series": 8, 
          "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
          "Aisc": 0.00279, "Aimp": -0.0003, "Bvoco": -0.01608, "Mbvoc": 0, "Bvmpo": -0.01608, 
          "Mbvmp": 0, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
          "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.0287, "B1": -0.0108 , 
          "B2": 0.001, "B3": -4E-05 , "B4": 5E-07, "B5": -3E-09, "DTC": 3, "FD": 1, "A": -3.47, 
          "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
          "C4": 1, "C5": 0, "C6": 1, "C7": 0} # B parameters assessed via regression of I AM modifier values
# =============================================================================
# #Alberino P
# module = {"Vintage": 2020, "Material": "mc-Si", "Cells_in_Series": 8, 
#           "Parallel_Strings": 2, "Isco": 3.5, "Voco": 5.36, "Impo": 3.3, "Vmpo": 4.568, 
#           "Aisc": 0.00279, "Aimp": -0.0003, "Bvoco": -0.01608, "Mbvoc": 0, "Bvmpo": -0.01608, 
#           "Mbvmp": 0, "N": 1, "IXO": 3.5, "IXXO": 2.05, "A0": 0.9645, "A1": 0.02753, 
#           "A2": -0.002848, "A3": -0.0001439, "A4": 0.00002219, "B0": 1.069, "B1": -0.0251, 
#           "B2": 0.0022, "B3": -8E-05, "B4": 1E-06, "B5": -6E-09, "DTC": 3, "FD": 1, 
#           "A": -3.47, "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
#           "C4": 1, "C5": 0, "C6": 1, "C7": 0}

class Photovoltaic():
         def __init__(self, latitude, longitude, altitude, timezone, time, dni, ghi, dhi, temp_amb, wind_amb, pressure):
            
            self.temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['close_mount_glass_glass']
                           
            self.times = time
            
            self.clearsky_dni = dni
            self.clearsky_ghi = ghi
            self.clearsky_dhi = dhi
                
            self.solpos = pvlib.solarposition.get_solarposition(self.times, latitude, longitude)
            
            self.dni_extra = pvlib.irradiance.get_extra_radiation(self.times)        
            
            self.airmass = pvlib.atmosphere.get_relative_airmass(self.solpos['apparent_zenith'])
                    
            self.am_abs = pvlib.atmosphere.get_absolute_airmass(self.airmass, pressure)
            
            self.aoi = pvlib.irradiance.aoi(39, 180,
                                   self.solpos['apparent_zenith'], self.solpos['azimuth'])
    
            self.total_irrad = pvlib.irradiance.get_total_irradiance(39, 180,
                                                            self.solpos['apparent_zenith'],
                                                            self.solpos['azimuth'],
                                                            self.clearsky_dni, self.clearsky_ghi, self.clearsky_dhi,
                                                            dni_extra=self.dni_extra,
                                                            model='haydavies')
            
            
            self.tcell = pvlib.temperature.sapm_cell(self.total_irrad['poa_global'],
                                            temp_amb, wind_amb,
                                            -2.98,-0.0471,1)
            
            self.effective_irradiance = pvlib.pvsystem.sapm_effective_irradiance(
            self.total_irrad['poa_direct'], self.total_irrad['poa_diffuse'],
            self.am_abs, self.aoi, module)
    
            self.dc = pvlib.pvsystem.sapm(self.effective_irradiance, self.tcell, module)
            self.annual_energy = self.dc['p_mp'].sum()
            
               
# column_values = ["Angle","Yield"]
# electrical_data = pd.DataFrame(data=energies, columns = column_values)
# electrical_data.fillna(0)
# electrical_data.loc['Total'] = electrical_data.select_dtypes(np.number).sum()
# pd.set_option('display.max_colwidth', 40)
# print(electrical_data)
# electrical_data.to_excel(r'D:\REM\SEM III\Project 2\03.06.2021\Results2.xlsx')



