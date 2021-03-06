#!/usr/bin/env python

"""
Stomatal optimisation based on xylem hydraulics (SOX) model.

SOX hypothesis: find the gs maximises the product of leaf photosynthesis
and the xylem hydraulic conductance, i.e. dA.K / dgs = 0

References
==========
* Eller CB, Rowland L, Oliveira RS, Bittencourt PRL, Barros FV, da Costa ACL,
  Meir P, Friend AD, Mencuccini M, Sitch S, et al. (2018). Modelling tropical
  forest responses to drought and El Niño with a stomatal optimization model
  based on xylem hydraulics. Philosophical transactions of the Royal Society
  of London. Series B, Biological sciences 373: 20170315.
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

def sox_optimisation(Vcmax25, Tleaf, Cs, PAR, press, psi_pd, p50, a_vuln,
                     rp_min, dq, LAI, k=0.5):

    C = CollatzC3()

    # Calculate Ci at the co-limitation point, i.e. the point where increasing
    # Ci won't lead to an increase in A because A would become light (Al) or
    # transport limited (Ae)
    Ci_col = C.calc_ci_at_colimitation_point(Cs, Tleaf, PAR, Vcmax25)

    # Calculate dCi
    dCi = Cs - Ci_col

    # Calculate dA
    A_Cs = C.calc_photosynthesis(Cs, Tleaf, PAR, Vcmax25) # Ambient leaf CO2
    A_col = C.calc_photosynthesis(Ci_col, Tleaf, PAR, Vcmax25) # co-limit point
    dA = A_Cs - A_col

    # Calculate dA/dCi
    dA_dci = dA / (dCi / press)

    # Calculate dK
    # The dK is computed as the linear gradient between K_pd and
    # K((psi_pd+p50)/2), see supplementary...this isn't in the main manuscript.
    # They use p50 as they argue this represents the steepest point of the
    # vulnerability func K(psi_m).
    K_pd = calc_xylem_hydraulic_conduc(psi_pd, p50, a_vuln)
    psi_mid = (psi_pd + p50) / 2.0
    K_psi_mid = calc_xylem_hydraulic_conduc(psi_mid, p50, a_vuln)
    dK = K_pd - K_psi_mid

    # Calculate dpsi_leaf
    dpsi_leaf = (psi_pd - psi_mid)

    # Calculate dK/dpsi_leaf
    dK_dpsi_leaf = dK / dpsi_leaf

    # plant hydraulic resistance (m2 s mol-1 H2O)
    rp = rp_min / K_pd

    # Calculate xi, the cost of stomatal opening in terms of loss of xylem
    # conductivity under low Ψpd and/or higher leaf-to-air vapour pressure.
    # Low xi = high hydraulic cost
    xi = 2.0 / ( (1.0 / K_pd) * dK_dpsi_leaf * rp * c.GSVGSC * dq)

    # Calculate gs at the co-limitation point
    gs_col = (A_col * press / (Cs - Ci_col)) * c.GSVGSC

    # Calculate gs (mol H2O m-2 s-1) using the analytical approximation for
    # the optimal SOX gs using the partial derivatives of A with respect to Ci
    # and K with respect to psi_m
    if dA_dci <= 0.0:
        gs = gs_col
    else:
        gs = 0.5 * dA_dci * (np.sqrt(4.0 * xi / dA_dci + 1) - 1.0)
        gs *= c.GSVGSC

    # Infer transpiration, assuming perfect coupling (mol H2O m-2 leaf s-1)
    E = dq * gs

    # Infer psi_leaf (MPa)
    psi_leaf = psi_pd - (rp_min / K_pd) * E

    # Infer plant hydraulic resistance (m2 s mol-1 H2O)
    K = calc_xylem_hydraulic_conduc( (psi_pd + psi_leaf) / 2.0, p50, a_vuln)
    rp = rp_min / K

    # Infer An (mol CO2 m-2 leaf s-1)
    gc = gs / c.GSVGSC
    An = C.calc_photosynthesis_given_gc(Cs, Tleaf, PAR, Vcmax25, gc, press)

    # Diagnose Ci (Pa)
    Ci = Cs - An / (gc / press)

    # Scale A (g C m-2 d-1) and E (mm d-1) up to the canopy using a
    # big-leaf approximation
    GPP = An * (1.0 - np.exp(-k * LAI)) / k
    GPP *= c.MOL_C_TO_GRAMS_C * c.SEC_TO_DAY
    E *= (1.0 - np.exp(-k * LAI)) / k
    E *= c.MOL_WATER_2_G_WATER * c.G_TO_KG * c.SEC_TO_DAY

    return (gs, GPP, E, Ci, rp)


def calc_xylem_hydraulic_conduc(psi, p50, a):
    """
    Calculate the normalised (0 to 1) xylem hydraulic conductance (K), computed
    from the inverse polynomial function that describes cavitation-induced
    embolism.

    Parameters:
    ----------
    psi : float
        water potential (MPa)
    p50 : float
        xylem pressure inducing 50% loss of hydraulic
        conductivity due to embolism (MPa)
    a : float
        shape of the vulnerability curve (unitless)

    Reference:
    ---------
    * Manzoni S et al. (2013) Hydraulic limits on maximum plant transpiration and
      the emergence of the safety – efficiency trade-off. New Phytol. 198,
      169-178.
    * Christoffersen BO et al. (2016) Linking hydraulic traits to tropical
      forest function in a size-structured and trait-driven model
      (TFS v.1-Hydro). Geosci. Model Dev. 9, 4227-4255.

    """
    return 1.0 / (1.0 + (psi / p50)**a)


if __name__ == "__main__":

    """
    Vcmax25 = 0.0001 # Maximum rate of rubisco activity 25C (mol m-2 s-1)
    Tleaf = 20.      # Leaf temp (deg C)
    Cs = 36.0        # leaf CO2 partial pressure (Pa)
    PAR = 0.002      # photosyn active radiation (mol m-2 s-1)
    press = 101325.0 # atmospheric pressure (Pa)
    psi_pd = -0.1    # pre-dawn soil water potential (MPa)
    p50 = -2.0       # xylem pressure inducing 50% loss of hydraulic
                     # conductivity due to embolism (MPa)
    a_vuln = 3.0     # shape of the vulnerability curve (unitless)
    rp_min = 2000.   # plant minimum hydraulic resistance
                     # (= 1/kmax; m2 s MPa mol-1 H2O)
    dq = 0.02474747  # Leaf to air vapor pressure deficit (unitless)
    LAI = 2.
    (gs, GPP, E, Ci, rp) = sox_optimisation(Vcmax25, Tleaf, Cs, PAR, press,
                                            psi_pd, p50, a_vuln, rp_min, dq,
                                            LAI)
    """

    Vcmax25 = 0.0001 # Maximum rate of rubisco activity 25C (mol m-2 s-1)
    Tleaf = 20.0 # Leaf temp (deg C)
    Cs = np.arange(1, 81) # leaf CO2 partial pressure (Pa)
    PAR = 0.002 # photosyn active radiation (mol m-2 s-1)
    press = 101325.0 # atmospheric pressure (Pa)
    psi_pd = -0.1    # pre-dawn soil water potential (MPa)
    p50 = -2.0       # xylem pressure inducing 50% loss of hydraulic
                     # conductivity due to embolism (MPa)
    a_vuln = 3.0     # shape of the vulnerability curve (unitless)
    rp_min = 2000.   # plant minimum hydraulic resistance
                     # (= 1/kmax; m2 s MPa mol-1 H2O)
    dq = 0.02474747  # Leaf to air vapor pressure deficit (unitless)
    LAI = 2.

    gs_store = []

    for i in range(80):

        (gs, GPP, E, Ci, rp) = sox_optimisation(Vcmax25, Tleaf, Cs[i], PAR,
                                                press, psi_pd, p50, a_vuln,
                                                rp_min, dq, LAI)
        gs_store.append(gs)

    Cs_umol_mol = Cs / press * c.MOL_TO_UMOL
    plt.plot(Cs_umol_mol, gs_store)
    plt.show()
