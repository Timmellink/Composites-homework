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
import th_effects as th
import epsilon_star_laminate as esl
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
alpha_array = asl.alpha_star_laminate(alph,layup)
# %% calculate Nth and Mth by using th_effects
NMth = th.th_effects(Cstar_lam,delta,alpha_array,z)
NMth = np.round(NMth, 7)
print("The thermal forces Nth and Mth are: "+str(NMth))
# %%  calculate epsilon and kappa
# first, calculate abd matrix
ABD_M = ab.ABD_matrix(Cstar_lam, z)
abd = np.linalg.inv(ABD_M)
NM= np.array([0,0,0,0,0,0]) # no load is applied
force = np.sum([NMth,NM],axis=0)
# then, calculate deformation, assuming there is no load applied
epsk = abd@force # material CS
epsk = np.round(epsk,4)
#%% calculate stress in each layer, in ply CS
# to do this, iterate over the layers, multiplying 
# [C*]*({eps*}-{alpha*}*deltaT)
# so, one needs the rotated deformation in ply CS
eps = epsk[0:3]
eps_star_lam = esl.epsilon_star_laminate(eps, layup)


