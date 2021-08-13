# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
from tespy.networks import Network
from tespy.components import (Sink, Source, Valve, Merge,
                              Subsystem, Compressor,
                              SolarCollector)
from tespy.connections import Connection

from tespy.tools.logger import logging
# from tespy.tools import logger
from tespy.tools.helpers import TESPyComponentError
import math

Temp_Units = {
    'C': [273.15, 1], 'F': [459.67, 5 / 9], 'K': [0, 1],
    'R': [0, 5 / 9]
}


class sdp_subsys(Subsystem):
    """
    subsystem definition for the sdp.

    It consists of the solar-roof tiles,
    connected in series and also includes a gap between the tiles, where air
    is drawn in. This model supposes, that the air is drawn from a ventilator
    at the end of the roof. The ventilator is not included in the subsystem.

    For connecting in series and parallel

    """

    def __init__(self, label,
                 num_sdp_series,
                 lkf_lin_sdp,
                 lkf_quad_sdp,
                 A_sdp,
                 Tamb,
                 E_sdp,
                 ks_sdp,
                 d_sdp,
                 l_sdp,
                 p_amb,
                 m_loss=None,
                 zeta_sdp=None):

        if not isinstance(label, str):
            msg = 'Subsystem label must be of type str!'
            logging.error(msg)
            raise ValueError(msg)

        elif len([x for x in [';', ', ', '.'] if x in label]) > 0:
            msg = 'Can\'t use ' + str([';', ', ', '.']) + ' in label.'
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.label = label

        if num_sdp_series <= 1:
            raise TESPyComponentError('Minimum number of sdp in series is 2.')
        else:
            self.num_sdp_series = num_sdp_series
        if m_loss is None and zeta_sdp is None:
            raise TESPyComponentError('Please provide either m_loss or zeta_sdp.')

        self.lkf_lin_sdp = lkf_lin_sdp
        self.lkf_quad_sdp = lkf_quad_sdp
        self.A_sdp = A_sdp
        self.Tamb = Tamb
        self.E_sdp = E_sdp
        self.ks_sdp = ks_sdp
        self.d_sdp = d_sdp
        self.l_sdp = l_sdp
        self.m_loss = m_loss
        self.zeta_sdp = zeta_sdp
        self.p_amb = p_amb

        self.comps = {}
        self.conns = {}
        self._create_comps()
        self._create_conns()

    def _create_comps(self):

        for i in range(self.num_sdp_series):
            j = str(i)

            self.comps['merge_' + j] = Merge(label=self.label +
                                             '_merge_' + j)

            self.comps['sdp_' + j] = SolarCollector(label=self.label +
                                                    '_sdp_' + j,
                                                    lkf_lin=self.lkf_lin_sdp,
                                                    lkf_quad=self.lkf_quad_sdp,
                                                    A=self.A_sdp,
                                                    E=self.E_sdp,
                                                    D=self.d_sdp,
                                                    L=self.l_sdp,
                                                    Tamb=self.Tamb,
                                                    ks=self.ks_sdp,
                                                    eta_opt=1)

            self.comps['valve_' + j] = Valve(label=self.label + '_Valve_' + j)

            if self.zeta_sdp is not None:
                self.comps['valve_' + j].set_attr(zeta=self.zeta_sdp)
            self.comps['source_' + j] = Source(label=self.label + '_Source_' + j)

    def _create_conns(self):

        for i in range(self.num_sdp_series):
            j = str(i)

            if i > 0:
                self.conns['mesd_' + j] = Connection(self.comps['merge_' +
                                                                str(i - 1)],
                                                     'out1',
                                                     self.comps['sdp_' + j],
                                                     'in1')

            self.conns['sdme_' + j] = Connection(self.comps['sdp_' + j],
                                                 'out1',
                                                 self.comps['merge_' + j],
                                                 'in1')
            
            # assign value for mass flow loss in each connection. 
            # This is done for modelling the Pressure Drop within the SRT system
            
            if self.m_loss is not None:
                if isinstance(self.m_loss, pd.DataFrame):                       # If multiple values for the leakage massflows are given, they are assigned accordingly for each Connection 
                   self.conns['sova_' + j] = Connection(
                        self.comps['source_' + j],
                        'out1',
                        self.comps['valve_' + j],
                        'in1',
                        p=self.p_amb,
                        T=self.Tamb,
                        fluid={'air': 1},
                        m=self.m_loss.iloc[i][1])
                else:
                    self.conns['sova_' + j] = Connection(                       # If there is a single value given for the leakage mass flow, it will be assigned for each Connection
                        self.comps['source_' + j],
                        'out1',
                        self.comps['valve_' + j],
                        'in1',
                        p=self.p_amb,
                        T=self.Tamb,
                        fluid={'air': 1},
                        m=self.m_loss)
            else:
                self.conns['sova_' + j] = Connection(                           # If m_loss is None m0 will be assigned for each connection 
                    self.comps['source_' + j],
                    'out1',
                    self.comps['valve_' + j],
                    'in1',
                    p=self.p_amb,
                    T=self.Tamb,
                    fluid={'air': 1},
                    m0=0.0001)                                                  # Starting value specification for mass flow 
                
            self.conns['vame_' + j] = Connection(
                self.comps['valve_' + j],
                'out1',
                self.comps['merge_' + j],
                'in2',
                p0=1.01)


class SDP_sucking():
    """
    Model of the SDP
    """

    def __init__(self,
                 sdp_in_series,
                 sdp_in_parallel,
                 loss_key_figures_lin=34,
                 loss_key_figures_quad=0.03,
                 width_sdp=0.244,
                 thickness_sdp=0.02,
                 l_sdp=0.3275,
                 ks=0.0001,
                 zeta=None):

        """
        Simulates the solar roof tile with a fan at the back of it, which is
        "sucking" air from the sdps and into the air-Source heat-pump.
        This TESPy model simulates each sdp in a series of sdps and thus
        includes effects of leakage and uneven temperature distribution.

        Parameters:
        -----------

            num_sdp_series : int
                number of SDPs connected in series for the model

            num_sdp_parallel : int
                number of SDPs connected in parallel for the model

            loss_key_figures_lin : float
                the losses of one sdp, linearly dependent on the temperature
                difference

            loss_key_figures_quad : float
                the losses of one sdp, quadratically dependent on the
                temperature difference

            width_sdp : float
                the width of the sdp in m

            thickness_sdp : float
                the thickness of the sdp in m

            l_sdp : float
                the length (or height) of the sdp in m
        """

        # -------------------------------------------------------------------
        # For the TESPy model...
        self.num_sdp_series = sdp_in_series
        self.num_sdp_parallel = sdp_in_parallel

        # loss figures of the sdp
        self.loss_key_figures_lin = loss_key_figures_lin
        self.loss_key_figures_quad = loss_key_figures_quad

        # width of the sdp in m
        self.width_sdp = width_sdp

        # thickness of the sdp in m
        self.thickness_sdp = thickness_sdp

        # summarized length of the sdp:
        self.l_sdp = l_sdp  # height/length of the sdp in m

        # equivalent diameter of the sdp to calculate the pressure loss via the
        # ks value
        self.d_sdp = math.sqrt(thickness_sdp * width_sdp * 2 / math.pi)

        # the surface area of the sdp in m²
        self.surface_area_sdp = width_sdp * l_sdp

        # zeta value for the leakage in each sdp
        self.zeta = zeta

        # ks value for the pressure loss inside each sdp
        self.ks = ks

    # %% initialisation of the sdp

    def init_sdp(self,
                 ambient_temp,
                 absorption_incl,
                 inlet_temp,
                 mass_flow,
                 p_amb=1.01325,
                 m_loss=None,
                 zeta=None,
                 print_res=True):
        """
        Initializing the TESPy simulation of the SDP. This initialization is
        stored in the folder "sdp", which is then used by the calculate_sdp
        function.

        Parameters:
        -----------

            ambient_temp : float
                ambient temperature to calculate the losses of the sdp

            inlet_temp : float
                temperature of the air at the inlet of the first sdp

            absorption_incl : float
                effective global irradiation on the sdp (optical losses already
                included!)

            mass_flow : float
                mass flow of air at the point behind the fan. It will evenly be
                devided by the model to each parallel string of sdps

            zeta : float
                the zeta value of the gap between the sdps to simulate the
                "leakage" of air, or for this model the stream of air drawn in
                between the sdps.

            m_loss : float
                the massflow of air at each "leakage" to initialize the zeta
                value. This is to be taken with care and does not really work
                yet, because the massflow cannot be the same for each sdp in
                the series of sdps. It is an optional parameter - better not
                use it for now.

            print_res : boolean, default = True
                if set True, it will print out the results of the simulation
                Set this parameter to False to suppress this.

        Returns:
        --------

            t_out : float
                The temperature at the outlet, behind the fan.

            p_fan : float
                The power of the fan needed to provide the massflow.

            m_out : float
                The massflow at the outlet, behind the fan. It should be the
                same, as what was given as mass_flow.
        """

        # network
        fluid_list = ['air']
        self.nw = Network(fluids=fluid_list,
                          p_unit='bar',
                          T_unit='C',
                          v_unit='m3 / s',
                          m_range=[0.0001, 11],
                          p_range=[0.9, 1.02],
                          h_range=[0, 100000],
                          T_range=[-10, 100])

        # %% components

        # Sinks & Sources
        feed_coll = Source('to collector')
        from_coll = Sink('from collector')

        # fan
        fan = Compressor('fan')

        mass_flow = mass_flow / self.num_sdp_parallel

        # sdp_subsystem
        sdp_sub = sdp_subsys('sdp',
                             num_sdp_series=self.num_sdp_series,
                             lkf_lin_sdp=self.loss_key_figures_lin,
                             lkf_quad_sdp=self.loss_key_figures_quad,
                             A_sdp=self.surface_area_sdp,
                             Tamb=ambient_temp,
                             E_sdp=absorption_incl,
                             ks_sdp=self.ks,
                             d_sdp=self.d_sdp,
                             l_sdp=self.l_sdp,
                             zeta_sdp=zeta,
                             m_loss=m_loss,
                             p_amb=p_amb)

        # %% Connections

        out_num = self.num_sdp_series - 1
        feedc_sdp = Connection(feed_coll,
                               'out1',
                               sdp_sub.comps["sdp_0"],
                               'in1')
        sdp_fan = Connection(sdp_sub.comps["merge_" + str(out_num)],
                             'out1',
                             fan,
                             'in1')
        self.nw.add_conns(feedc_sdp, sdp_fan)

        fan_fromc = Connection(fan, 'out1', from_coll, 'in1')
        self.nw.add_conns(fan_fromc)

        self.nw.add_subsys(sdp_sub)

        fan.char_warnings = False
        fan.set_attr(
            eta_s=0.5,
            design=['eta_s'],
            offdesign=['eta_s_char']
        )

        # %% Connection parameters
        feedc_sdp.set_attr(
            m0=mass_flow,
            T=inlet_temp,
            fluid={'air': 1},
            p=p_amb,
        )

        sdp_fan.set_attr(
            p0=1,
        )

        fan_fromc.set_attr(
            p=p_amb,
            m=mass_flow,
        )

        # %% solving

        mode = 'design'
        self.nw.check_network()
        # self.nw.save('sdp')
        self.nw.solve(mode=mode)
        self.nw.save('sdp')
        if print_res:
            self.nw.print_results()

        # %% return characteristic parameters

        t_out = fan_fromc.T.val_SI - Temp_Units[self.nw.T_unit][0]
        p_fan = fan.P.val_SI * self.num_sdp_parallel
        m_out = fan_fromc.m.val_SI * self.num_sdp_parallel
        return t_out, p_fan, m_out

    # %% calculation of the sdp

    def calculate_sdp(self,
                      ambient_temp,
                      absorption_incl,
                      inlet_temp,
                      mass_flow,
                      ks_SRT,
                      p_amb=1.01325,
                      print_res=True):
        """
        TESPy calculation of the sdp after initialization.

        Parameters:
        -----------

            ambient_temp : float
                ambient temperature to calculate the losses of the sdp

            inlet_temp : float
                temperature of the air at the inlet of the first sdp

            absorption_incl : float
                effective global irradiation on the sdp (optical losses already
                included!)

            mass_flow : float
                mass flow of air at the point behind the fan. It will evenly be
                devided by the model to each parallel string of sdps
                
            ks_SRT: float
                ks/roughness value for SRTs, delivered to calculate the pressure 
                Drop over a hole String of SRTs in Off-Design mode.
                The goal of using the ks value is to adjust the pressure drop of 
                the SRTs in an iterative way so that they fit the more accure 
                results for the pressure drop from the CFD modelling.

        Returns:
        --------

            t_out : float
                The temperature at the outlet, behind the fan.

            p_fan : float
                The power of the fan needed to provide the massflow.

            m_out : float
                The massflow at the outlet, behind the fan. It should be the
                same, as what was given as mass_flow.
        """

        if ambient_temp is not None:
            for comp in self.nw.comps['object']:
                if isinstance(comp, SolarCollector):
                    comp.set_attr(Tamb=ambient_temp)

        if absorption_incl is not None:
            for comp in self.nw.comps['object']:
                if isinstance(comp, SolarCollector):
                    comp.set_attr(E=absorption_incl)
                    
        if ks_SRT is not None:
            for comp in self.nw.comps['object']:
                if isinstance(comp, SolarCollector):
                    comp.set_attr(ks=ks_SRT)

        if inlet_temp is not None:

            for conn in self.nw.conns['object']:
                a = type(conn)
                if isinstance(conn.source, Source):
                    conn.set_attr(T=inlet_temp)

        if mass_flow is not None:
            mass_flow = mass_flow / self.num_sdp_parallel
            for conn in self.nw.conns['object']:
                if isinstance(conn.source, Compressor):
                    conn.set_attr(m=mass_flow)

        self.nw.save("sdp")
        # %% solving
        mode = 'offdesign'
        self.nw.solve(mode=mode, design_path='sdp')
        if print_res:
            self.nw.print_results()

        # get parameters
        for conn in self.nw.conns['object']:
            if isinstance(conn.source, Compressor) & isinstance(conn.target, Sink):
                t_out = conn.T.val_SI - Temp_Units[self.nw.T_unit][0]
                m_out = conn.m.val_SI * self.num_sdp_parallel

        for comp in self.nw.comps['object']:
            if isinstance(comp, Compressor):
                p_fan = comp.P.val * self.num_sdp_parallel
                
        return t_out, p_fan, m_out
    
    
    def plot_temperature_curve(self, 
                               p_amb):
        
        """
        Calculates & plots the pressure Drop within the SRT String 
        
        """
        
        Valve_name = "Valve {}"                                                 # Creating variable for iterating through Valves 1 - 12 
        i = 0
        p_Valve = []
        for comp in self.nw.comps['object']:
            if isinstance(comp, Valve):
                
                if i == 0:
                    p_ref_Valve0 = comp.pr.val                                  # The pressure within valve 1 is taken as the reference pressure for the pressure drop calculation within the SRT String
                    print('\nPressure in Valve 1 (reference pressure of the system) in Bar:', comp.pr.val.round(8),     # Extracting pressure value from components bin within the network. Extracted is the Pressure in Bar for Valve 1
                          '\nCumulative (!!) pressure loss for the SRT string in Pa:\n')
                                    
                p_Valve.append(                                                 # Append the list for technical values for the Valves
                {
                    'Valve': Valve_name.format(str(i+1)),                       # Valve label  and pressure difference [Pa]
                    'Pressure [bar]': comp.pr.val,                              # Pressure within Valve [bar]
                    'Pressure Difference [Pa]': ((p_ref_Valve0 - comp.pr.val)*-1)*100000    #Pressure Difference [Pa] of Valve x compared to the pressure in Valve 1. Negative because of a loss of pressure in the system. 
                }
                )
                if i > 0:
                    print('Pr from', Valve_name.format(str(i)), 'to', Valve_name.format(str(i+1)), 'in [Pa]:', (((p_ref_Valve0 - comp.pr.val)*-1)*100000).round(2)) #Print pressure Drop from each Valve to the next Valve.
                i=i+1
                
        P_Valve=pd.DataFrame(p_Valve)                                           # Save results in DataFrame
        
        fig, ax = plt.subplots(figsize=(8,5), sharex=True)                      # Plotting results

        ax.plot(P_Valve['Valve'], P_Valve['Pressure [bar]'], linestyle='solid', color='red', label='Pressure [bar]')
        ax.set_xlabel("Valve", fontsize=14)
        ax.set_ylabel("Pressure [bar]", fontsize=14)
        ax2= ax.twinx()
        ax2.plot(P_Valve['Valve'], P_Valve['Pressure Difference [Pa]'], linestyle='solid', label='Pressure Difference [Pa]')
        ax2.set_ylabel("Druckdifferenz [Pa]", fontsize=14)
        plt.title('Pressure Drop & Pressure Difference to ambient Pressure in the Valves of the subsystem')
        fig.legend(bbox_to_anchor=[0.875, 0.85], loc="upper right")
        plt.gcf().autofmt_xdate() 
         
        
        return P_Valve
 