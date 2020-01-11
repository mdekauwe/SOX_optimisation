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
    * Clark DB, Mercado LM, Sitch S, Jones CD, Gedney N, Best MJ, Pryor M,
      Rooney GG, Essery RLH, Blyth E, et al. 2011. The Joint UK Land
      Environment Simulator (JULES), Model description – Part 2: Carbon fluxes
      and vegetation. Geoscientific Model Development Discussions 4: 641–688.

    """
    # Note in Clark et al. Ko25=30.0*1E4, using 1E3 to match Eller, check if
    # that is a mistake
    def __init__(self, Oi=210.0, gamstar25=42.75, Kc25=30.0, Ko25=30.0*1E3,
                 Ec=79430.0, Eo=36380.0, theta_hyperbol=0.9995,
                 quantum_yield=0.3, absorptance=0.8, gs_model=None, gamma=0.0,
                 Eav=60000.0, deltaSj=650.0, deltaSv=650.0, Hdv=200000.0,
                 Q10_Kc=2.1, Q10_Ko=1.2, Q10_Vcmax=2.0, Tlower=10.0,
                 Tupper=50.0):


        self.Oi = Oi
        self.Kc25 = Kc25 # Pa
        self.Ko25 = Ko25 # Pa
        self.Eav = Eav
        self.deltaSv = deltaSv
        self.Hdv = Hdv
        self.Q10_Ko = Q10_Ko # Q10 value for MM constants for O2
        self.Q10_Kc = Q10_Kc # Q10 value for MM constants for CO2
        self.Q10_Vcmax = Q10_Vcmax # Q10 value for carboxylation of Rubisco
        self.Tlower = Tlower # Lower temperature for carboxylation
        self.Tupper = Tupper # Upper temperature for carboxylation


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
        Vcmax = self.correct_vcmax_for_temperature(Vcmax25, Tleaf)

        # Michaelis Menten constants for CO2 (Pa)
        Kc = self.Q10_func(self.Kc25, self.Q10_Kc, Tleaf)

        # Michaelis Menten constants for O2 (Pa)
        Ko = self.Q10_func(self.Ko25, self.Q10_Ko, Tleaf)

        print(Vcmax, Kc, Ko)


    def Q10_func(self, k25, Q10, Tleaf):
        """
        Q10 function to calculate parameter change with temperature
        """

        return k25 * (Q10**((Tleaf - 25.0) / 10.0))

    def correct_vcmax_for_temperature(self, Vcmax25, Tleaf):
        """
        Correct Vcmax based on defined by PFT-specific upper and lower
        temperature params, see Clark et al. (mol CO2 m-2 s-1)
        """
        num = self.Q10_func(Vcmax25, self.Q10_Vcmax, Tleaf)
        den = (1.0 + math.exp(0.3 * (Tleaf - self.Tupper))) * \
                (1.0 + math.exp(0.3 * (self.Tlower - Tleaf)))

        return num / den


if __name__ == "__main__":

    # Maximum rate of rubisco activity 25C (mol m-2 s-1)
    Vcmax25 = 0.0001 #
    Tleaf = 35.0

    Cs = 400.0

    # photosynthetically active radiation (mol m-2 s-1)
    PAR = 0.002

    C = CollatzC3()
    C.calc_photosynthesis(Cs, Tleaf, PAR, Vcmax25)
