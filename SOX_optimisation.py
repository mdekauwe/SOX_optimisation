#!/usr/bin/env python

"""
Stomatal optimisation based on xylem hydraulics (SOX) model.

Reference
=========
* Eller, C.B., Rowland, L., Mencuccini, M., Rosas, T., Williams, K., Harper, A.,
  Medlyn, B.E., Wagner, Y., Klein, T., Teodoro, G.S., Oliveira, R.S.,
  Matos, I.S., Rosado, B.H., Fuchs, K., Wohlfahrt, G., Montagnani, L., Meir, P.,
  Sitch, S. and Cox, P.M. (2020), Stomatal optimisation based on xylem
  hydraulics (SOX) improves land surface model simulation of vegetation
  responses to climate. New Phytol., DOI: doi:10.1111/nph.16419.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (11.01.2020)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from collatz_photosynthesis import CollatzC3
import constants as c

def main(Vcmax25, Tleaf, Cs, PAR, press, psi_pd, p50, a_vuln, rp_min, dq):

    Ci = Cs * 0.7

    C = CollatzC3()

    Ci_col = C.calc_ci_at_colimitation_point(Ci, Tleaf, PAR, Vcmax25)

    # Calculate dA/dCi
    Ci1 = Cs
    Ci2 = Ci_col
    dCi = Ci1 - Ci2

    A1 = C.calc_photosynthesis(Ci1, Tleaf, PAR, Vcmax25)
    A2 = C.calc_photosynthesis(Ci2, Tleaf, PAR, Vcmax25)
    dA = A1 - A2

    dA_dci = dA / (dCi / press)

    # Calculate dV/dLWP numerically
    V1 = cavitation_func(psi_pd, p50, a_vuln)
    psi_mid = (psi_pd + p50) / 2.0
    V2 = cavitation_func(psi_mid, p50, a_vuln)
    dV = V1 - V2

    dV_dLWP_V = (dV / (psi_pd - psi_mid)) * (1.0 / V1)

    # Calculate xi
    rp_d = c.GSVGSC * (rp_min / V1) * dq / 2.0
    xi = 1.0 / (rp_d * dV_dLWP_V)

    print(xi)

def cavitation_func(P, P50, a):
    return 1.0 / (1.0 + (P / P50)**a)


if __name__ == "__main__":

    #N = 100
    #Vcmax25 = 0.0001 # Maximum rate of rubisco activity 25C (mol m-2 s-1)
    #Tleaf = np.repeat(20, 100) # Leaf temp (deg C)
    #Ca = 36                    # leaf CO2 partial pressure (Pa)
    #PAR = np.repeat(0.002, 100) # photosyn active radiation (mol m-2 s-1)
    #press = np.repeat(101325., 100) #  atmospheric pressure (Pa)

    #N  100
    Vcmax25 = 0.0001 # Maximum rate of rubisco activity 25C (mol m-2 s-1)
    Tleaf = 20.      # Leaf temp (deg C)
    Cs = 36.0        # leaf CO2 partial pressure (Pa)
    PAR = 0.002      # photosyn active radiation (mol m-2 s-1)
    press = 101325.0 # atmospheric pressure (Pa)
    psi_pd = -0.1    # pre-dawn soil water potential (MPa)
    p50 = -2.0       # xylem pressure inducing 50% loss of hydraulic
                     # conductivity due to embolism (MPa)
    a_vuln = 3.0     # shape of the vulnerability curve (unitless)
    rp_min = 2000.   # minimum plant hydraulic resistance
                     # (= 1/kmax; mol-1 m s MPa)
    dq = 0.02474747  # Leaf to air vapor pressure deficit (unitless)
    main(Vcmax25, Tleaf, Cs, PAR, press, psi_pd, p50, a_vuln, rp_min, dq)
