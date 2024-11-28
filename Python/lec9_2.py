# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:25:05 2024

@author: Timme
"""
# %% import packages
import numpy as np
import ply_edges as pl
import Cstar_laminate as cs
import stiffness as s
import ABD as ab
import alpha_star as ap
# %% set parameters
E1 = 128e9
E2 = 9e9
G12 = 5e9
nu12 = 0.35
t = 0.15e-3 # thickness ply
alph1 = 2e-6
alph2 = 50e-6
alph = np.array([alph1, alph2, 0]) # vector of thermal coefficients in material CS
delta = 100
layup = [0,45,90,-45]
n = len(layup)*2
# %% Calculate alpha star cell array
# first calculate C, then calculate Cstar laminate
# Then, one can calculate alpha star cell array
C = s.stiffness(E1,E2,nu12,G12)
Cstar_lam = cs.Cstar_laminate(layup,C,n)
z = pl.ply_edges(t,n)
alpha_array = ap.alpha_star(alph,layup)


# %%
