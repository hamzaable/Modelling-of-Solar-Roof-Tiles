# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 22:06:44 2020

@author: Vaishnavi Phadke
"""
from tespy.networks import network
from tespy.components import (sink, source, valve, merge,
                              subsystem, compressor,
                              solar_collector)
from tespy.connections import connection

from tespy.tools.logger import logging
# from tespy.tools import logger
from tespy.tools.helpers import TESPyComponentError
import math

class sdp_subsys(subsystem):
    """
    subsystem definition for the sdp.

    It consists of the solar-roof tiles,
    connected in series and also includes a gap between the tiles, where air
    is drawn in. This model supposes, that the air is drawn from a ventilator
    at the end of the roof. The ventilator is not included in the subsystem.

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
            self.comps['merge_' + j] = merge(label=self.label +
                                             '_merge_' + j)
            self.comps['sdp_' + j] = solar_collector(label=self.label +
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
            self.comps['valve_' + j] = valve(label=self.label + '_valve_' + j)
            if self.zeta_sdp is not None:
                self.comps['valve_' + j].set_attr(zeta=self.zeta_sdp)
            self.comps['source' + j] = source(label=self.label + '_source_' + j)

    def _create_conns(self):

        for i in range(self.num_sdp_series):
            j = str(i)

            if i > 0:
                self.conns['mesd_' + j] = connection(self.comps['merge_' +
                                                                str(i - 1)],
                                                     'out1',
                                                     self.comps['sdp_' + j],
                                                     'in1')

            self.conns['sdme_' + j] = connection(self.comps['sdp_' + j],
                                                 'out1',
                                                 self.comps['merge_' + j],
                                                 'in1')
            if self.m_loss is not None:
                self.conns['sova_' + j] = connection(
                    self.comps['source' + j],
                    'out1',
                    self.comps['valve_' + j],
                    'in1',
                    p=self.p_amb,
                    T=self.Tamb,
                    fluid={'air': 1},
                    m=self.m_loss)
            else:
                self.conns['sova_' + j] = connection(
                    self.comps['source' + j],
                    'out1',
                    self.comps['valve_' + j],
                    'in1',
                    p=self.p_amb,
                    T=self.Tamb,
                    fluid={'air': 1},
                    m0=0.0001)
            self.conns['vame_' + j] = connection(
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
        "sucking" air from the sdps and into the air-source heat-pump.
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
        self.d_sdp = math.sqrt(thickness_sdp*width_sdp*2/math.pi)

        # the surface area of the sdp in m²
        self.surface_area_sdp = width_sdp*l_sdp

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
        self.nw = network(fluids=fluid_list,
                          p_unit='bar',
                          T_unit='C',
                          v_unit='m3 / s',
                          m_range=[0.0001, 11],
                          p_range=[0.9, 1.02],
                          h_range=[0, 100000],
                          T_range=[-10, 100])

        # %% components

        # sinks & sources
        feed_coll = source('to collector')
        from_coll = sink('from collector')

        # fan
        fan = compressor('fan')

        mass_flow = mass_flow/self.num_sdp_parallel

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

        # %% connections

        out_num = self.num_sdp_series-1
        feedc_sdp = connection(feed_coll,
                               'out1',
                               sdp_sub.comps["sdp_0"],
                               'in1')
        sdp_fan = connection(sdp_sub.comps["merge_"+str(out_num)],
                             'out1',
                             fan,
                             'in1')
        self.nw.add_conns(feedc_sdp, sdp_fan)

        fan_fromc = connection(fan, 'out1', from_coll, 'in1')
        self.nw.add_conns(fan_fromc)

        self.nw.add_subsys(sdp_sub)

        # %% component parameters
        fan.set_attr(
                eta_s=0.5,
                design=['eta_s'],
                offdesign=['eta_s_char']
                )

        # %% connection parameters
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
        t_out = fan_fromc.T.val_SI-self.nw.T[self.nw.T_unit][0]
        p_fan = fan.P.val_SI*self.num_sdp_parallel
        m_out = fan_fromc.m.val_SI*self.num_sdp_parallel
        return t_out, p_fan, m_out

    # %% calculation of the sdp

    def calculate_sdp(self,
                      ambient_temp,
                      absorption_incl,
                      inlet_temp,
                      mass_flow,
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
            for comp in self.nw.comps.index:
                if isinstance(comp, solar_collector):
                    comp.set_attr(Tamb=ambient_temp)

        if absorption_incl is not None:
            for comp in self.nw.comps.index:
                if isinstance(comp, solar_collector):
                    comp.set_attr(E=absorption_incl)

        if inlet_temp is not None:
            for conn in self.nw.conns.index:
                if isinstance(conn.source, source):
                    conn.set_attr(T=inlet_temp)

        if mass_flow is not None:
            mass_flow = mass_flow/self.num_sdp_parallel
            for conn in self.nw.conns.index:
                if isinstance(conn.source, compressor):
                    conn.set_attr(m=mass_flow)

        self.nw.save("sdp")
        # %% solving
        mode = 'offdesign'
        self.nw.solve(mode=mode, design_path='sdp')
        if print_res:
            self.nw.print_results()

        # get parameters
        for conn in self.nw.conns.index:
            if isinstance(conn.source, compressor) & isinstance(conn.target, sink):
                t_out = conn.T.val_SI-self.nw.T[self.nw.T_unit][0]
                m_out = conn.m.val_SI*self.num_sdp_parallel

        for comp in self.nw.comps.index:
            if isinstance(comp, compressor):
                p_fan = comp.P.val*self.num_sdp_parallel

        return t_out, p_fan, m_out