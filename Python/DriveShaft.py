# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:03:08 2025

@author: Timme
"""

# %% import statements
import CompositeProperties as cp
import numpy as np

# %% set parameters
E1 = 120e9
E2 = 8e9
G12 = 3.5e9
nu12 = 0.35
alph1 = 0.16e-6
alph2=27.6e-6

# %% Calculate Nxy and G
C = cp.stiffness(E1,E2,nu12,G12) # calculate stiffness of single ply
Cst = cp.# calculate rotated stifness matrix
# Calculate ABD
# Calculate abd
# Gxy = 1/(a66*h)
# calculate Nxy = 2T/piD^2
# calculate epsxy0 = a66*Nxy
# calculate theta = L/R*epsxy0