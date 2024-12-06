# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:14:37 2024

@author: Timme
"""
# %% import statements
import composite_properties as cp
import numpy as np
import thermal_effects as te
import importlib as il

# %% Tensile test variables
L = 250e-3 # m
W = 25e-3 # m
layup = [45,0,-45,90]*3
n = len(layup)*2
E1 = 130e9 # Pa
E2 = 10e9
nu12 = 0.33
G12 = 5.1e9 
t = 0.15e-3 
# strength values
S1t = 2250e6
S1c = 1300e6
S2t = 110e6
S2c = 130e6
S6 = 205e6
# %% calculate abd
il.reload(cp)
C = cp.stiffness(E1,E2,nu12,G12)# calculate C
Cstar_laminate = cp.Cstar_laminate(layup,C)# Calculate C star laminate
z = cp.ply_edges(t,n)# Calculate ply edges z
ABD = cp.ABD_matrix(Cstar_laminate,z)# calculate ABD
abd = np.linalg.inv(ABD)# calcualte abd
# %% print abd
for i in abd:
    print("\n")
    for j in i:
        print(f"{j:01} ", end="")