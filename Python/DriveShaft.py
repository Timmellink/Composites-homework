# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:03:08 2025

@author: Timme
"""

# %% import statements
import CompositeProperties as cp
import numpy as np
import TensileTest as tt
import FailTest as ft
from importlib import reload

# %% set parameters
E1 = 120e9
E2 = 8e9
G12 = 3.5e9
nu12 = 0.35
alph1 = 0.16e-6
alph2=27.6e-6
layup = [45,-45]
t=0.1e-3 # test thickness
T = 500 # max torque
R = 50e-3 # radius shaft
D = 2*R 
L = 400e-3
theta = np.radians(0.5)
S1t = 1400e6
S1c = 1200e6
S2t = 80e6
S2c = 110e6
S6 = 100e6
strength = [S1c, S1t, S2c, S2t, S6]

# %% Calculate Nxy and G
C = cp.stiffness(E1,E2,nu12,G12) # calculate stiffness of single ply
Cst = cp.Cstar_laminate(layup,C) # calculate rotated stifness matrix
n = len(layup)*2

h = n*t
#Gxy = 1/(abd[5][5]*h) # Gxy = 1/(a66*h)
Nxy = 2*T/(np.pi*D**2) # calculate Nxy = 2T/piD^2
epsxy0 = theta*R/L # calculate epsxy0 = theta * R/L

# %% calculate t
t = Nxy/(4*Cst[1][2][2]*epsxy0) # t = A66/(4*C166*)
print("The minimum needed thickness is "+str(t))

# %% Calculate ABD to see if there is shear extension coupling
z = cp.ply_edges(t,n)
ABD = cp.ABD_matrix(Cst,z) # Calculate ABD
abd = np.linalg.inv(ABD) # Calculate abdK

# %% Calculate stresses
#reload(tt)
reload(ft)
NM = [0, 0, Nxy, 0, 0, 0] # calculate NM
# epUKsk0 = abd@NM # Calculate epsk0
# sig_star = tt.CalculateStress(NM, Cst, z)# calculate sig_star
# sig_mat = tt.RotateMaterial(sig_star,layup)# rotate to sig_mat
fails = ft.TsaiTest(NM,layup, strength,Cst, t)
# %%
