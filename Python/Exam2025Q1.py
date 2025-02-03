# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 08:42:03 2025

@author: Timme
"""

# %% import modules
import numpy as np
import CompositeProperties as cp
import thermal_effects as te
from importlib import reload
import FailTest as ft

# %% set composite properties

To = 190
Tf = 20
delta = Tf-To
t = 0.2e-3
layup = [0,0,90,90]

# material properties

E1 = 100e9
E2 = 10e9
G12 = 5e9
nu12 = 0.3
alph1 = 0.2e-6
alph2 = 50e-6

# strength data

S1t = 2100e6
S1c = 1300e6
S2t = 100e6
S2c = 120e6
S6 = 200e6

# %% calculate C and ABD

C = cp.stiffness(E1, E2, nu12, G12)
Cstar_lam = cp.Cstar_laminate(layup, C)
n = len(layup)*2
z = cp.ply_edges(t, n)
ABD = cp.ABD_matrix(Cstar_lam, z)

# %% calculate fictive thermal forces and defs.

alpha_vec = [alph1, alph2, 0]
alpha_star_lam = te.alpha_star_laminate(alpha_vec, layup)
Fth = te.calculate_NM_thermal(Cstar_lam, delta, alpha_star_lam, z)
abd = cp.inv(ABD)
eps_k = abd@Fth
print(f"the strains and curvatures are : {eps_k}")

# %% plot material stresses
sig_th = te.ThermalStress(Fth, alpha_star_lam, delta, Cstar_lam, z)
sig_mat = cp.RotateMaterial(sig_th, layup)
te.PlotThermalStress(sig_mat, 0, z)
te.PlotThermalStress(sig_th, 0, z)

# %% max. stress criterion
strength = [S1c, S1t, S2c, S2t, S6]
fails = ft.MaxStressTestPlyFail(sig_th, strength)

# %% claculate Ex and Gxy for sample A
a11 = abd[0][0]
a66 = abd[2][2]
h = t*n
Ex_A = 1/(a11*h)
Gxy_A = 1/(a66*h)
print(f"The value for Ex_A : {Ex_A*1e-9} GPa, while Gxy_A : {Gxy_A*1e-9} GPa")
# %% calculate Ex_B and Gxy_B
layup2 = [angle-45 for angle in layup]
C_B = cp.stiffness(E1, E2, nu12, G12)
C_star_lam_B = cp.Cstar_laminate(layup2, C_B)
ABD_B = cp.ABD_matrix(C_star_lam_B, z)
abd_B = cp.inv(ABD_B)
a11_B = abd_B[0][0]
a66_B = abd_B[2][2]

Ex_B = 1/(a11_B*h)
Gxy_B = 1/(a66_B*h)
print(f"the value for Ex_B : {Ex_B*1e-9} GPa, while Gxy_B : {Gxy_B*1e-9} GPa")