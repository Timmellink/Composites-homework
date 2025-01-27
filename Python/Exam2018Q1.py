# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:14:51 2025

@author: Timme
"""

# %% import statements

import numpy as np
import CompositeProperties as cp
import thermal_effects as te
from importlib import reload
# %% set Composite properties

E1 = 100e9
E2 = 10e9
G12 = 3.5e9
nu12 = 0.3
layup = [45]*4+[-45]*4
t = 0.15e-3
alph1 = 1e-6
alph2 = 50e-6
deltaT = -100
# %% calculate abd

reload(cp)
C = cp.stiffness(E1, E2, nu12, G12)
Cstar_lam = cp.Cstar_laminate(layup, C)
n = cp.nPlies(layup)
z = cp.ply_edges(t, n)
ABD = cp.ABD_matrix(layup, z)
abd = cp.inv(ABD)
# %% calculate NMth and eps
alpha_vec = [alph1, alph2, 0]
alpha_star = te.alpha_star_laminate(alpha_vec, layup)
NMth = te.calculate_NM_thermal(Cstar_lam, deltaT, alpha_star, z)
