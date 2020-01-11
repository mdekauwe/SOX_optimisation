#!/usr/bin/env python

"""
Collatz photosynthesis model.


That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (11.01.2020)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import math
import constants as c

class CollatzC3(object):

    """
    Collatz photosynthesis model.


    Reference
    =========
    * Collatz, G. J., Ball, J. T., Grivet, C., and Berry, J. A. (1991)
      Physiological and environmental regulation of stomatal conductance,
      photosynthesis and transpiration: amodel that includes alaminar boundary
      layer, Agr. Forest Meteorol., 54, 107–136.

    """
    # Note in Clark et al. Ko25=30.0*1E4, using 1E3 to match Eller, check if
    # that is a mistake
    def __init__(self, Oi=210.0, gamstar25=42.75, Kc25=30.0, Ko25=30.0*1E3,
                 Ec=79430.0, Eo=36380.0, theta_hyperbol=0.9995,
                 quantum_yield=0.3, absorptance=0.8, gs_model=None, gamma=0.0,
                 Eav=60000.0, deltaSj=650.0, deltaSv=650.0, Hdv=200000.0,
                 Q10_Kc=2.1, Q10_Ko=1.2, Q10=2.0, Tlower=10.0, Tupper=50.0):


        self.Oi = Oi
        self.Kc25 = Kc25 # Pa
        self.Ko25 = Ko25 # Pa
        self.Eav = Eav
        self.deltaSv = deltaSv
        self.Hdv = Hdv
        self.Q10_Ko = Q10_Ko
        self.Q10_Kc = Q10_Kc
        self.Q10 = Q10
        self.Tlower = Tlower
        self.Tupper = Tupper


    def calc_photosynthesis(self, Cs=None, Tleaf=None, Par=None, Vcmax25=None):
        """
        Parameters
        ----------
        Cs : float
            leaf surface CO2 concentration (mol mol-1)
        Tleaf : float
            leaf temp (deg C)
        PAR : float
            photosynthetically active radiation (mol m-2 s-1)
        Vcmax25 : float
            Maximum rate of rubisco activity 25C (mol m-2 s-1)

        """
        Tk = Tleaf + c.DEG_2_KELVIN

        # Max rate of rubisco activity (mol m-2 s-1)
        Vcmax = (Vcmax25 * (self.Q10**(0.1 * (Tleaf - 25.0))) /
                            ((1.0 + math.exp(0.3 * (Tleaf - self.Tupper))) *
                            (1.0 + math.exp(0.3 * (self.Tlower - Tleaf)))))

        # Michaelis Menten constants for CO2 (Pa)
        Kc = self.Q10_func(self.Kc25, self.Q10_Kc, Tleaf)

        # Michaelis Menten constants for O2 (Pa)
        Ko = self.Q10_func(self.Ko25, self.Q10_Ko, Tleaf)

        print(Vcmax, Kc, Ko)


    def Q10_func(self, k25, Q10, Tleaf):
        """ Q10 function to calculate parameter change with temperature """

        return k25 * (Q10**((Tleaf - 25.0) / 10.0))


    def peaked_arrh(self, k25, Ea, Tk, deltaS, Hd):
        """ Temperature dependancy approximated by peaked Arrhenius eqn,
        accounting for the rate of inhibition at higher temperatures.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC or 298 K
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]
        deltaS : float
            entropy factor [J mol-1 K-1)
        Hd : float
            describes rate of decrease about the optimum temp [J mol-1]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.

        """
        arg1 = self.arrh(k25, Ea, Tk)
        arg2 = 1.0 + np.exp((298.15 * deltaS - Hd) / (298.15 * c.RGAS))
        arg3 = 1.0 + np.exp((Tk * deltaS - Hd) / (Tk * c.RGAS))

        return arg1 * arg2 / arg3

    def arrh(self, k25, Ea, Tk):
        """ Temperature dependence of kinetic parameters is described by an
        Arrhenius function.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC or 298 K
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.
        """
        return k25 * np.exp((Ea * (Tk - 298.15)) / (298.15 * c.RGAS * Tk))


if __name__ == "__main__":

    # Maximum rate of rubisco activity 25C (mol m-2 s-1)
    Vcmax25 = 0.0001 #
    Tleaf = 35.0

    Cs = 400.0

    # photosynthetically active radiation (mol m-2 s-1)
    PAR = 0.002

    C = CollatzC3()
    C.calc_photosynthesis(Cs, Tleaf, PAR, Vcmax25)
