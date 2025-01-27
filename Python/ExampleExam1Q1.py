# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 08:44:44 2025

@author: Timme
"""

#%% import statements

import numpy as np
import CompositeProperties as cp
import thermal_effects as te
import FailTest as ft
from importlib import reload
# %% set properties

layup = [0] + [45]*2+[90]*2+[-45]*2+[0]
h = 0.15e-3
Tcure = 200
Troom = 20
delta = Troom-Tcure
E1 = 100e9
E2 = 10e9
G12 = 5e9
nu12 = 0.3
alph1 = 0.2e-6
alph2 = 50e-6
S1t = 2100e6
S1c = 1400e6
S2t = 100e6
S2c = 120e6
S6 = 200e6
# %% calculate abd matrix

reload(cp)
n = len(layup)*2 # calculate n
t = h/n # calculate t 
C = cp.stiffness(E1, E2, nu12, G12)
Cstar_lam = cp.Cstar_laminate(layup, C)
z = cp.ply_edges(t, n)
ABD = cp.ABD_matrix(Cstar_lam, z)
abd = np.linalg.inv(ABD)
# %% calculate fictive thermal forces

reload(te)
alpha = [alph1, alph2, 0] # calculate alpha vector
alpha_star_lam = te.alpha_star_laminate(alpha, layup)# calculate alpha_star lam
NMth = te.calculate_NM_thermal(Cstar_lam, delta, alpha_star_lam, z)
# %% calculate deformations
def_therm = NMth@abd # multiply NMth with abd
# %% calculate thermal stresses

reload(ft)
strength = [S1c, S1t, S2c, S2t, S6] # calculate strength
reload(te)
sig_therm = te.ThermalStress(NMth, alpha_star_lam, delta, Cstar_lam, z)
sig_therm_mat = cp.RotateMaterial(sig_therm, layup)
te.PlotThermalStress(sig_therm, 0, z) # plot stresses in ply CS
fails = ft.MaxStressTestPlyFail(sig_therm_mat, strength)# test material stresses against max stress
# %% calculate E modulus and G modulus

a11 = abd[0][0]
a66 = abd[2][2]
Ex = 1/(a11*h) # Ex = 1/(a11*h)
Gxy = 1/(a66*h) # Gxy = 1/(a66*h)
print(f"The E modulus in x direction is : {Ex*1e-9} GPa, while the shear modulus is : {Gxy*1e-9} GPa.")