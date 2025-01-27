# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:14:37 2024

@author: Timme
"""
# %% import statements
import CompositeProperties as cp
import numpy as np
import thermal_effects as te
import importlib as il
import matplotlib.pyplot as plt
import FailTest as ft
from importlib import reload

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
h=t*n

# strength values
S1t = 2280e6
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

# %% calculate values engineering constants
Ex = 1/(abd[0][0]*h)# Ex = 1/(a11*h) h is laminate thickness
Ey = 1/(abd[1][1]*h)# Ey = 1/(a22*h)
Gxy = 1/(abd[2][2]*h)# Gxy = 1/(a33*h)
nuxy = -(abd[1][0])/(abd[0][0]) # nuxy = -a21/a11

print("Ex : ",str(Ex))
print("Ey : ",str(Ey))
print("Gxy : ",str(Gxy))
print("nu_xy : ",str(nuxy))

# %% plot ply CS stresses
Fx = 50e3 # Fx = 50kN
Nx = Fx/W # determine NM vector
NM = [Nx, 0, 0, 0, 0, 0]
cp.PlotPlyStress(NM, Cstar_laminate,z,0)

# %% Test whether plies fail
strength = (S1c,S1t,S2c,S2t,S6)
Fx = 50e3 # Fx = 50kN
Nx = Fx/W # determine NM vector
NM = [Nx, 0, 0, 0, 0, 0]
fail = ft.TsaiTest(NM,layup,strength, Cstar_laminate, t)
n_lst = np.array(range(n))
#Failures,MaterialStresses = RotateMaterial(stresses,layup)
#print("The plies that fail are number(s) #", n_lst[Failures==True])
stress= cp.CalculateStress(NM,Cstar_laminate,z)
MatStr = cp.RotateMaterial(stress,layup)
print("The plies that fail are number(s) #", n_lst[fail==True])

# %% plot material stresses
strength = (S1c,S1t,S2c,S2t,S6)
Fx = 51.9e3 # Fx = 51.9kN
Nx = Fx/W # determine NM vector
NM = [Nx, 0, 0, 0, 0, 0]
cp.PlotMatStr(NM,Cstar_laminate,z,0,layup)
cp.PlotMatStr(NM,Cstar_laminate,z,1,layup)
cp.PlotMatStr(NM,Cstar_laminate,z,2,layup)
