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

def main(Vcmax25, Tleaf, Cs, PAR):

    Ci = Cs * 0.7

    C = CollatzC3()
    An = C.calc_photosynthesis(Ci, Tleaf, PAR, Vcmax25)

    print(An)

if __name__ == "__main__":

    Vcmax25 = 0.0001 # Maximum rate of rubisco activity 25C (mol m-2 s-1)
    Tleaf = 35.0     # Leaf temp (deg C)
    Cs = 40.0       # leaf CO2 partial pressure (Pa)
    PAR = 0.002      # photosynthetically active radiation (mol m-2 s-1)

    main(Vcmax25, Tleaf, Cs, PAR)
