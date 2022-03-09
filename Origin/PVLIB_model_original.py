# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 22:06:11 2020

@author: Vaishnavi Phadke
"""
import pvlib
import numpy as np
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
          "B2": 0.0010, "B3": - 3.5338E-05 , "B4": 5.3488E-07, "B5": -2.9287E-09, "DTC": 3, "FD": 1, "A": -3.47, 
          "B": -0.0594, "C0": 1, "C1": 0, "C2": 0.912848952156834, "C3": 0.0582212364987667, 
          "C4": 1, "C5": 0, "C6": 1, "C7": 0, "celltype": "monoSi", "gamma_pmp": -0.3792} # B parameters assessed via regression of I AM modifier values

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
         def __init__(self, latitude, longitude, altitude, timezone, time, dni, ghi, dhi, temp_amb, wind_amb, pressure, module_number):
            
            self.temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['close_mount_glass_glass']
                           
            self.times = time
            
            self.dni = dni
            self.clearsky_dni = dni
            self.clearsky_ghi = ghi
            self.clearsky_dhi = dhi
            
            self.module_number = module_number
                
            self.solpos = pvlib.solarposition.get_solarposition(self.times, latitude, longitude)
            
            self.dni_extra = pvlib.irradiance.get_extra_radiation(self.times)        
            
            self.airmass = pvlib.atmosphere.get_relative_airmass(self.solpos['apparent_zenith'])
                    
            self.am_abs = pvlib.atmosphere.get_absolute_airmass(self.airmass, pressure)
            
            self.aoi = pvlib.irradiance.aoi(39, 171,
                                   self.solpos['apparent_zenith'], self.solpos['azimuth'])
            
            #Manipulation of poa_direct(PVGIS) (E_dir on horizontal plane) to DNI-tilted for get_total_irradiance calculation
            self.dni = self.dni/np.cos(np.radians(90-self.solpos['apparent_elevation'])) 
            
            self.total_irrad = pvlib.irradiance.get_total_irradiance(39, 171,
                                                            self.solpos['apparent_zenith'],
                                                            self.solpos['azimuth'],
                                                            self.dni, self.clearsky_ghi, self.clearsky_dhi,
                                                            dni_extra=self.dni_extra,
                                                            model='haydavies')
            
            
            # Estimating cell temperatue via Faimann Model
            # Faiman, D. (2008). “Assessing the outdoor operating temperature of 
            # photovoltaic modules.” Progress in Photovoltaics 16(4): 307-315.
            self.tcell = pvlib.temperature.faiman(self.total_irrad['poa_global'],
                                                        temp_amb,
                                                        wind_amb,
                                                        u0=17.896,               # Combined heat loss factor coefficient [W*m^-2*C^-1]
                                                        u1=2.015)                # Combined heat loss factor influenced by wind [W*m^-2*C^-1(m/s)]
            
            self.effective_irradiance = pvlib.pvsystem.sapm_effective_irradiance(
            self.total_irrad['poa_direct'], self.total_irrad['poa_diffuse'],
            self.am_abs, self.aoi, module)
            
    
            #self.dc = pvlib.pvsystem.sapm(self.effective_irradiance, self.tcell, module)
            #self.annual_energy = self.dc['p_mp'].sum()
            
            #Estimates parameters for the CEC single diode model (SDM) using the SAM SDK
            #For Details on the single module specs see module definition in SDP.py 
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
                self.tcell,                                                     # cell temperature 
                alpha_sc=module['Aisc'],                                   # See SDP.py
                a_ref=self.mp_fit_cec_sam[4],                                   # The product of the usual diode ideality factor n (unitless), number of cells in series, and cell thermal voltage at reference conditions [V]
                I_L_ref=self.mp_fit_cec_sam[0],                                 # The light-generated current (or photocurrent) at reference conditions [A]
                I_o_ref=self.mp_fit_cec_sam[1],                                 # The dark or diode reverse saturation current at reference conditions [A]
                R_sh_ref=self.mp_fit_cec_sam[3],                                # The shunt resistance at reference conditions, in ohms.
                R_s=self.mp_fit_cec_sam[2],                                     # The series resistance at reference conditions, in ohms.
                EgRef=1.121,                                                    # The energy bandgap at reference temperature in units of eV. 1.121 eV for crystalline silicon. EgRef must be >0. For parameters from the SAM CEC module database, EgRef=1.121 is implicit for all cell types in the parameter estimation algorithm used by NREL
                dEgdT=- 0.0002677,                                              # The temperature dependence of the energy bandgap at reference conditions in units of 1/K. May be either a scalar value (e.g. -0.0002677 as in 1) or a DataFrame (this may be useful if dEgdT is a modeled as a function of temperature). For parameters from the SAM CEC module database, dEgdT=-0.0002677 is implicit for all cell types in the parameter estimation algorithm used by NREL.
                irrad_ref=1000,                                                 # Reference irradiance in W/m^2
                temp_ref=25)                                                    # Reference cell temperature in C
                        
            
            #Solving the single-diode equation to obtain a photovoltaic IV curve. Important is the output p_mp
            self.single_diod_IV = pvlib.pvsystem.singlediode(
                photocurrent=self.single_diod_parameters[0],                    # Light-generated current IL (photocurrent)
                saturation_current=self.single_diod_parameters[1],              # Series resistance Rs under desired IV curve conditions
                resistance_series=self.single_diod_parameters[2],               # Series resistance Rs under desired IV curve conditions
                resistance_shunt=self.single_diod_parameters[3],                # Shunt resistance Rsh under desired IV curve conditions.
                nNsVth=self.single_diod_parameters[4],                          # The product of three components: 1) the usual diode ideality factor n, 2) the number of cells in series Ns, and 3) the cell thermal voltage Vth. 
                ivcurve_pnts=None,                                              # Number of points in the desired IV curve. If None or 0, no points on the IV curves will be produced.
                method='lambertw')                                              # Determines the method used to calculate points on the IV curve. The options are 'lambertw', 'newton', or 'brentq'
            
            # Storing module Power to pass it on to the Excel Sheet 
            # Power output calculation for mutliple modules (DC)
            self.annual_energy = self.module_number*self.single_diod_IV['p_mp'].sum() #in Watt hours, because p_mp is given in Watt
            
            # Power output calculation for mutliple modules (AC) with PVWatts Inverter Method (NREL)
            """
                
                NREL’s PVWatts inverter model.
                The PVWatts inverter model [1] calculates inverter efficiency η as a function of input DC power
                [1] A. P. Dobos, “PVWatts Version 5 Manual,” http://pvwatts.nrel.gov/downloads/pvwattsv5.pdf (2014).
                
                For the estimation of the fitting inverter see Project2_further-calculations.xlsx Sheet 'Inverter DC/AC'
                Rated power of System = 12 * 16 * 15 W = 2.880 kW --> pdc0 assumed to be 3.5 kW
            
                Parameters:
                -----------
                
                    pdc : numeric
                        DC power. Same unit as pdc0.
                        
                    pdc0 : numeric
                        DC input limit of the inverter. Same unit as pdc.
                        
                    eta_inv_nom : numeric
                        Nominal inverter efficiency. [unitless]. Default 0.96
                        
                    eta_inv_ref : numeric
                        Reference inverter efficiency. PVWatts defines it to be 0.9637 and is included here for flexibility. [unitless]
                        The reference inverter effiency_inv_ref is taken from the CEC data for the actual most typical inverter, which is 0.9637 (Pac0/Pdc0)
                    
                   
                Outputs:
                --------
                
                    power_ac : numeric
                        AC power. Same unit as pdc0.
                                            
            """         
            
            self.power_ac = pvlib.inverter.pvwatts(
                pdc=self.annual_energy,                                         
                pdc0= 6000,                                                      # Used Inverter is SMA Sunny Boy 1.5            
                eta_inv_nom=0.961,                                               # Maximum efficiency 97,2 % / Euro-eta 96,1 % 
                eta_inv_ref=0.9637)                                              # based on CEC Data
                
                   



