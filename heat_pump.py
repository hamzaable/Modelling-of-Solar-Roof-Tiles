#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:18:14 2019

@author: christian
"""
import CoolProp.CoolProp as cp

class HeatPump():
    """
    Heat pump object.
    """

    def __init__(self,
                 p_el_max,
                 v_max=3500.,
                 deltaT_max=5.,
                 hp_type="analytic_07",
                 min_uptime=60.,
                 ):
        """
        Heat pump object.

        Parameters:
        -----------

            p_el_max : float
                maximum electrical power provided by the heat pump

            t_out : float
                Minimum outlet temperature of the air-stream, leaving the heat
                pump

            v_max : float
                Maximaler Volumenstrom der Wärmepumpe

            hp_type : string
                Type of the heat-pump

        """

        self.p_el_max = p_el_max
        self.deltaT_max = deltaT_max
        self.v_max = v_max
        self.hp_type = hp_type
        self.min_uptime = min_uptime
        self.time_left = 0.
        self.is_on = False
        self.hp_type = hp_type
        # prepare the link to the correct cop-function
        if self.hp_type == "analytic_07":
            self.x = 482.9
            self.y = -1.39
            self.calc_cop = self.calc_cop_analytic
        elif self.hp_type == "analytic_10" or self.hp_type == "analytic_14":
            self.x = 2222.31
            self.y = -1.83
            self.calc_cop = self.calc_cop_analytic
        elif self.hp_type == "simple_model":
            self.calc_cop = self.calc_cop_sm

    def refresh_timeleft(self, t):
        """
        Refresh the timeleft from min_uptime.
        
        By either starting from beginning, or subtracting t from time_left.
        """
        if not self.is_on:
            self.time_left = self.min_uptime-t
            self.is_on = True
        else:
            self.time_left = self.time_left - t
            if self.time_left < 0: self.time_left = 0
            self.is_on = True

    def calc_cop_staffel(self,
                    t_source,
                    t_cons):
        """
        simple heat pump model for air-source heat-pumps derived from the paper
        "A review of domestic heat pumps" by Iain Staffel et. al,
        Energy Environ. Sci., 2012. DOI: 10.1039/c2ee22653g

        Parameters:
        -----------

            t_cons : float
                Temperature required by the consumer in °C

            t_source : float
                Temperature given in from the source used by the heat pump

        Returns:
        --------

            cop : float
                the COP of the heat pump
        """

        delta_t = t_cons - t_source
        cop = (6.81 - 0.121 * (delta_t) + 0.00063 * (delta_t)**2)
        return cop

    def calc_cop_ruhnau(self,
                    t_source,
                    t_cons,
                    hp_type="ashp"):
        """
        heat pump regression model "Ruhnau".

        Simple heat pump regression curves for the COP of
        ashp: air source heat pumps,
        gshp: ground source heat pumps,
        wshp: ground-water source heat pumps,
        derived from the paper
        O. Ruhnau, L. Hirth, und A. Praktiknjo,
        „Time series of heat demand and heat pump efficiency for energy system
        modeling“,
        Sci Data, Bd. 6, Nr. 1, S. 189, Dez. 2019
        doi: 10.1038/s41597-019-0199-y.

        Parameters:
        -----------

            t_cons : float
                Temperature required by the consumer in °C

            t_source : float
                Temperature given in from the source used by the heat pump

            hp_typ : string, default: ashp
                Type of heat pump to use (ashp, gshp, wshp)

        Returns:
        --------

            cop : float
                the COP of the heat pump
        """

        delta_t = t_cons - t_source
        if hp_type == "ashp":
            cop = (6.08 - 0.09 * (delta_t) + 0.0005 * (delta_t)**2)
        elif hp_type == "gshp":
            cop = (10.29 - 0.21 * (delta_t) + 0.0012 * (delta_t)**2)
        elif hp_type == "wshp":
            cop = (9.97 - 0.2 * (delta_t) + 0.0012 * (delta_t)**2)
        return cop


    def calc_cop_analytic(
            self,
            t_source,
            t_usage):
        """
        simple heat pump model for air-source heat-pumps derived from data from
        Datasheets by the enterprise Wolf.

        Parameters:
        -----------

            t_usage : float
                Temperature of the consumer in °C

            t_source : float
                Outlet temperature from the sdp

        Returns:
        --------

            cop : float
                the cop of the heat pump
        """

        delta_t = t_usage - t_source
        cop = self.x * delta_t ** self.y
        return cop

    def available_heat_of_input(
            self,
            t_sdp,
            mass_flow):
        """
        calculates the available thermal energy from the incoming air.

        The available energy is calculated with the assumption, that the
        heat pump is able to cool down and use the energy of a certain delta T
        in the air with a certain mass flow - regardless of the temperature.
        It neglects any effects like icing etc.
        """
        c_p = cp.PropsSI('C', 'P', 101325, 'T', t_sdp+273.15, 'Air')
        q_therm = mass_flow*c_p*(self.deltaT_max)
        return q_therm

    def q_th(
            self,
            t_in,
            t_usage,
            p_hp_el=None
            ):
        """
        Calculate the available thermal power output of the heat pump.

        at a given temperature difference.

        Parameters
        ----------
            t_in: float
                temperature used by the heat pump

            t_usage: float
                desired temperature for the end use case.

        Returns
        -------
            q_th: float
                energy output of the heat pump.

        """
        if p_hp_el is None:
            p_hp_el = self.p_el_max
        cop = self.calc_cop(t_in, t_usage)
        q_th = p_hp_el * cop
        return q_th
