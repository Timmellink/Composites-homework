# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 08:43:43 2025

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
h = 2e-3
layup = [45, -45]
D = 60e-3
L = 100e-2

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

# %% calculate ABD
C = cp.stiffness(E1, E2, nu12, G12)
Cstar_lam = cp.Cstar_laminate(layup, C)
n = len(layup)*2
t = h/n
z = cp.ply_edges(t, n)
ABD = cp.ABD_matrix(Cstar_lam, z)
abd = cp.inv(ABD)
# %% inbetween steps
A = ABD[0:3,0:3]
a = abd[0:3,0:3]
Re = np.array([[1,0,0],[0,1,0],[0,0,2]])
Re_inv = cp.inv(Re)
T = cp.transformation(layup[0])
Tinv = cp.inv(T)
K = Re@T@Re_inv
K_inv = cp.inv(K)

# %% determine N vector
eps1 = 0.002
eps = [eps1, 0, 0]
eps_star = K_inv@eps
N = A@eps_star
Nxy = N[2]

# %% determine torque
To = np.pi*D**2/2*Nxy# T = piD^2/2*A61*Kinv*eps1
print(f"the torque is : {To} Nm")

# %% max stress criterion
To_max = 8000
Nxy_max = 2*To_max/(np.pi*D**2)
NM_max = [0, 0, Nxy_max, 0, 0, 0]
epsk_max = abd@NM_max
sig_max_star = cp.CalculateStress(NM_max, Cstar_lam, z)
sig_max_mat = cp.RotateMaterial(sig_max_star, layup)
strength = [S1c,S1t,S2c,S2t,S6]
fails = ft.MaxStressTest(NM_max, layup, strength, Cstar_lam, t)