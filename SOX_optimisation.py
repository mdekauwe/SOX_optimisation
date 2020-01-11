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

def main(Vcmax25, Tleaf, Cs, PAR, press):

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

    print(dA_dci)

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

    main(Vcmax25, Tleaf, Cs, PAR, press)
