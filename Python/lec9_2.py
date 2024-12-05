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
import alpha_star_laminate as asl
import thermal_effects as th
import epsilon_star_laminate as esl
import matplotlib.pyplot as plt

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
h = t*n
x = np.linspace(0, h, n + 1)  # height of laminate subdivided in n+1 elements

# %% Calculate alpha star cell array
# first calculate C, then calculate Cstar laminate
# Then, one can calculate alpha star cell array
C = s.stiffness(E1,E2,nu12,G12)
Cstar_lam = cs.Cstar_laminate(layup,C,n)
z = pl.ply_edges(t,n)
alpha_array = asl.alpha_star_laminate(alph,layup)

# %% calculate Nth and Mth by using th_effects
NMth = th.thermal_effects(Cstar_lam,delta,alpha_array,z)
NMth = np.round(NMth, 7)
print("The thermal forces Nth and Mth are: "+str(NMth))

# %%  calculate epsilon and kappa
# first, calculate abd matrix
ABD_M = ab.ABD_matrix(Cstar_lam, z)
abd = np.linalg.inv(ABD_M)
NM= np.array([0,0,0,0,0,0]) # no load is applied
force = np.sum((NMth,NM),axis=0)
# then, calculate deformation, assuming there is no load applied
epsk0 = abd@force # material CS
epsk0 = np.round(epsk0,6)
print("The deformations are: "+str(epsk0))

#%% calculate stress in each layer, in ply CS
# to do this, iterate over the layers, multiplying 
# [C*]*({eps*}-{alpha*}*deltaT)
# so, one needs the rotated deformation in ply CS
eps0 = epsk0[0:3]
kappa = epsk0[3:6]
stresses = []
for i in range(n):
    sig_k_b = Cstar_lam[i]@(eps0+z[i]*kappa-alpha_array[i]*delta)
    sig_k_e = Cstar_lam[i]@(eps0+z[i+1]*kappa-alpha_array[i]*delta)
    stresses.append([sig_k_b,sig_k_e])
stresses = np.array(stresses)

# %% plot
str1 = stresses[:,:,0] # stresses in 1 direction
str1_lie = np.reshape(str1, (-1))  # get the stresses in one row
xr = x[::-1] # reverse the x-coordinates
xd = [[i,i]for i in xr] # repeat the elements twice
xd_lie = np.reshape(xd,-1) # reshape to one row
plt.plot(str1_lie,xd_lie[1:-1]) #plot from top to bottom (only take begin and end once)
