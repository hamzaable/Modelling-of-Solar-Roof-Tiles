# -*- coding: utf-8 -*-

import pvlib


# from pvlib.pvsystem import pvsystem
# pip install NREL-PySAM
# NREL PySAM has to be installed to run fit_cec_sam


class Photovoltaic():
    def __init__(self, latitude, longitude, altitude, timezone, m_azimut, m_tilt, module_number, time, dni, ghi, dhi, albedo, a_r, irrad_model, module, temp_amb, wind_amb, pressure,cell_temp ):
    
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
        
        self.albedo = albedo
        self.a_r = a_r
        self.irrad_model = irrad_model
        self.module = module
            
        #Radiation, power & module temperature calculations
        
        # Calculates solar position according to location and time 
        self.solpos = pvlib.solarposition.get_solarposition(
            self.times,                                                     # Must be localized or UTC will be assumed
            self.latitude,                                                  # Latitude in decimal degrees. Positive north of equator, negative to south
            self.longitude)                                                 # Longitude in decimal degrees. Positive east of prime meridian, negative to west
        
        # Determine extraterrestrial radiation from day of year
        self.dni_extra = pvlib.irradiance.get_extra_radiation(
            self.times)                                                     # Timestamp
        
        # Calculate relative (not pressure-adjusted) airmass at sea level
        self.airmass = pvlib.atmosphere.get_relative_airmass(
            self.solpos['apparent_zenith'])
                
        # Determine absolute (pressure-adjusted) airmass from relative 
        # airmass and pressure.
        self.am_abs = pvlib.atmosphere.get_absolute_airmass(
            self.airmass,                                                   # relative airmass
            pressure)                                                       # ambient pressure                           
        
        self.aoi = pvlib.irradiance.aoi(
            self.m_tilt,                                                    # Panel tilt from horizontal
            self.m_azimut,                                                  # Panel azimuth from north
            self.solpos['apparent_zenith'], 
            self.solpos['azimuth'])
        
        # Determine total in-plane irradiance and its beam, sky diffuse and ground reflected components, 
        # using the specified sky diffuse irradiance model.
        self.total_irrad = pvlib.irradiance.get_total_irradiance(
            self.m_tilt,                                                    
            self.m_azimut,                                                  
            self.solpos['apparent_zenith'],
            self.solpos['azimuth'],
            self.clearsky_dni,                                              # Direct Normal Irradiance
            self.clearsky_ghi,                                              # Global horizontal irradiance
            self.clearsky_dhi,                                              # Diffuse horizontal irradiance
            dni_extra=self.dni_extra,                                       # extraterrestrial radiation
            model=self.irrad_model,                                         # Irradiance model, in this case haydavies
            albedo=self.albedo)                                             # see SDP.py
    
                    
        #Estimating cell temperatue via Faimann Model
        self.tcell = cell_temp                                                    # Combined heat loss factor influenced by wind [(W/m^2)/(C)]
        
        #Estimating Angle of incidence modifier(IAM) using the Martin and Ruiz for diffuse radiation
        self.IAM_mod_diff = pvlib.iam.martin_ruiz_diffuse(
            self.m_tilt,                                                    # Surface tilt angles in decimal degrees. The tilt angle is defined as degrees from horizontal (e.g. surface facing up = 0, surface facing horizon = 90) surface_tilt must be in the range [0, 180]
            a_r=self.a_r,                                                   # The angular losses coefficient. This is an empirical dimensionless parameter. Values of a_r are generally on the order of 0.08 to 0.25 for flat-plate PV modules. a_r must be greater than zero.
            c1=0.4244,                                                      # First fitting parameter for the expressions that approximate the integral of diffuse irradiance coming from different directions. c1 is given as the constant 4 / 3 / pi (0.4244) 
            c2=None)                                                        # Second fitting parameter for the expressions that approximate the integral of diffuse irradiance coming from different directions. If c2 is None, it will be calculated according to the linear relationship in IEC 61853-3 
        
        #Calculating SAPM spectral loss coefficient, F1.
        self.F1 = pvlib.pvsystem.sapm_spectral_loss(
            self.am_abs,                                                    # Absolute airmass
            self.module)                                                    # A dict, Series, or DataFrame defining the SAPM performance parameters. Defined in SDP.py
                    
        #Estimating Angle of incidence modifier(IAM) using the Martin and Ruiz for direct radiation
        self.IAM_mod_dir = pvlib.iam.martin_ruiz(
            self.aoi,                                                       # The angle of incidence between the module normal vector and the sun-beam vector in degrees.
            a_r=self.a_r)                                                   # The angular losses coefficient. This is an empirical dimensionless parameter. Values of a_r are generally on the order of 0.08 to 0.25 for flat-plate PV modules. a_r must be greater than zero.
        
        self.effective_irradiance  = self.F1 * (self.total_irrad['poa_direct'] * self.IAM_mod_dir 
                                                + self.IAM_mod_diff[0] * self.total_irrad['poa_sky_diffuse']
                                                + self.IAM_mod_diff[1] * self.total_irrad['poa_ground_diffuse'])   
                    
        """
            
            Manual Effective Irradiances Effective Irradiance calculation 
            --> Ground reflection considered
    
            Ee = f_1(AM_a) (E_b f_2(AOI) + f_d E_d)"
            or
            Ee = F1 * (poa_direct * F2 + module['FD'] * poa_diffuse)
        
            Parameters:
            -----------
            
                F1 : numeric
                    The SAPM spectral loss coefficient
                
                poa_direct : OrderedDict or DataFrame
                    Direct irradiance on module plane 
                    
                poa_sky_diffuse : OrderedDict or DataFrame
                    sky Diffuse irradiance on module plane 
                    
                poa_ground_diffuse : OrderedDict or DataFrame
                    ground Diffuse irradiance on module plane 
                    
                self.IAM_mod_dir : numeric
                    Angle of incidence modifier(IAM) for direct radiation
                    
                self.IAM_mod_diff[0] : numeric
                    incidence angle modifiers (iam) for diffuse sky irradiance
                    
                self.IAM_mod_diff[1] : numeric
                    incidence angle modifiers (iam) for ground-reflected irradiance
                
                                        
        """         
        
        #Estimates parameters for the CEC single diode model (SDM) using the SAM SDK
        #For Details on the single module specs see module definition in SDP.py 
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
        
                            
        #Calculating parameters for solving the single diod model by De Soto et. al model
        self.single_diod_parameters = pvlib.pvsystem.calcparams_desoto(
            self.effective_irradiance, 
            self.tcell,                                                     # cell temperature 
            alpha_sc=self.module['Aisc'],                                   # See SDP.py
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
        # Power output calculation for mutliple modules
        self.annual_energy = self.module_number*self.single_diod_IV['p_mp'].sum() #in Watt hours, because p_mp is given in Watt
    
        
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
        # Faiman, D. (2008). “Assessing the outdoor operating temperature of 
        # photovoltaic modules.” Progress in Photovoltaics 16(4): 307-315.
        self.tcell = pvlib.temperature.faiman(self.total_irrad['poa_global'],
                                                    temp_avg,
                                                    wind_amb,
                                                    u0=25.0,                # Combined heat loss factor coefficient [W*m^-2*C^-1]
                                                    u1=6.84)                # Combined heat loss factor influenced by wind [W*m^-2*C^-1(m/s)]
    
