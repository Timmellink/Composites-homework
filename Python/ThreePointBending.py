# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:19:13 2024

@author: Timme
"""

# %% import statements
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import CompositeProperties as cp

# %% import text file
widths = [6, 12, 14, 13]
df = pd.read_csv("3pb_data.txt", header=[0], sep='\s+')
#df2.columns = ["time", "disp", "force"]

# %% set variables
t = 0.15e-3 # ply thickness
n = 24 # number of plies
h = t*n # laminate thickness
Ls = 100e-3 # sample span length
b_s = 13e-3 # sample width
cd = {'E1' : 130e9,
      'E2' : 10e9,
      'G12' : 5.2e9,
      'nu12' : 0.33,
      'alp1' : 0.16e-6,
      'alp' : 27.6e-6
      } # composite data
sd = {'S1t' : 2280e6,
      'S1c' : 1300e6,
      'S2t' : 110e6,
      'S2c' : 130e6,
      'S6' : 205e6}# strength data 
layup1 = [45,0,-45,90]*3
layup2 = [0]*3+[45,90,-45,-45,90,45,45,90,-45]

# %% plot data
df.plot(x="disp",y="force")

# %% calculate slope m
df1 = df[df.disp==1] # find force at disp = 1.0
Fd1 = df1.force.squeeze() # squeeze() to get a scalar value
df3 = df[df.disp==3] # find force at disp = 3.0
Fd3 = df3.force.squeeze()
# calculate slope m = delta F/delta d
delta_F = Fd3-Fd1 # delta F = Fd3 - Fd1
delta_d = 3-1 # delta d = 3 - 1
m = delta_F/delta_d 


# %% calculate flexural stiffness
Efs = (Ls*1e-3)**3*m/(4*b_s*1e-3*(h*1e-3)**3) # Secant stiffness is L^3m/4bh^3 (see standard) 

# %% to compare flexural modulus, one has to calculate Efx
# This can be done by calculating abd
C = cp.stiffness(cd['E1'],cd['E2'],cd['nu12'],cd['G12']) # first, calculate the stiffness C of one ply 
C_star_laminate = cp.Cstar_laminate(layup1, C) # Then, calculate C* laminate, based on layup (45,0,-45,90)3s
z = cp.ply_edges(t, n) # calculate ply edges
ABD = cp.ABD_matrix(C_star_laminate, z)# Calculate ABD matrix
abd = np.linalg.inv(ABD) # calculate abd

# %%  Find Efx by analyzing d11
Efx = 12/(abd[3][3]*h**3) # Efx = 12/d11*h^3
Efy = 12/(abd[4][4]*h**3) # Efy = 12/d22*h^3

# %% repeat the procedure using second layup
Cstar_lam2 = cp.Cstar_laminate(layup2,C) # first, calculate C*laminate 2
ABD2 = cp.ABD_matrix(Cstar_lam2, z) # then calculate ABD using C*lam 2
abd2 = np.linalg.inv(abd) # Then, calculate abd
Efx2 = 12/(abd2[3][3]*h**3)
Efy2 = 12/(abd2[4][4]*h**3)

# %% compare the values
Efs_m = Efs*1e-6 # flexular stiffness in MPa
Efx_m, Efy_m =  (Efx*1e-6,Efy*1e-6) # bending stiffnesses layup1 in MPa
Efx2_m, Efy2_m = (Efx2*1e-6,Efy2*1e-6)