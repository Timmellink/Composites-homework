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
import FailTest as ft
#import TensileTest as tt
import thermal_effects as te
from importlib import reload

# %% import text file
df = pd.read_csv("3pb_data.txt", header=[0], sep='\s+')

# %% set variables
Ls = 100e-3 # sample span length (meters)
b_s = 13e-3 # sample width (meters)
cd = {'E1' : 130e9,
      'E2' : 10e9,
      'G12' : 5.2e9,
      'nu12' : 0.33,
      'alp1' : 0.16e-6,
      'alp2' : 27.6e-6,
      't' : 0.15e-3 # ply thickness
      } # composite data
# strength data
S1t =  2280e6
S1c = 1300e6
S2t = 110e6
S2c =  130e6
S6 =  205e6
strength = (S1c,S1t,S2c,S2t,S6)
layup1 = [45,0,-45,90]*3
layup2 = [0]*3+[45,90,-45,-45,90,45,45,90,-45]
n = len(layup1)*2 # number of plies
h = cd['t']*n # laminate thickness
AlphaVec = [cd['alp1'], cd['alp2'], 0] # alpha vector in material CS
delta = 23-143 # negative temperature change

# %% plot data
#df.plot(x="disp",y="force")

# %% calculate slope m
disp1 = 0.1; disp3 = 3.5
df1 = df[df.disp==disp1] # find force at disp = 1.0
Fd1 = df1.force.squeeze() # squeeze() to get a scalar value
df3 = df[df.disp==disp3] # find force at disp = 3.0
Fd3 = df3.force.squeeze()

# calculate slope m = delta F/delta d
delta_F = Fd3-Fd1 # delta F = Fd3 - Fd1
delta_d = disp3-disp1 # delta d = 3 - 1
m = delta_F/delta_d 


# %% calculate flexural stiffness
Efs = (Ls*1e3)**3*m/(4*b_s*1e3*(h*1e3)**3) # Secant stiffness is L^3m/4bh^3 (see standard) 

# %% to compare flexural modulus, one has to calculate Efx
# This can be done by calculating abd
C = cp.stiffness(cd['E1'],cd['E2'],cd['nu12'],cd['G12']) # first, calculate the stiffness C of one ply 
C_star_laminate = cp.Cstar_laminate(layup1, C) # Then, calculate C* laminate, based on layup (45,0,-45,90)3s
z = cp.ply_edges(cd['t'], n) # calculate ply edges
ABD = cp.ABD_matrix(C_star_laminate, z)# Calculate ABD matrix
abd = np.linalg.inv(ABD) # calculate abd

# Find Efx by analyzing d11
Efx = 12/(abd[3][3]*h**3) # Efx = 12/d11*h^3
Efy = 12/(abd[4][4]*h**3) # Efy = 12/d22*h^3

# repeat the procedure using second layup
Cstar_lam2 = cp.Cstar_laminate(layup2,C) # first, calculate C*laminate 2
ABD2 = cp.ABD_matrix(Cstar_lam2, z) # then calculate ABD using C*lam 2
abd2 = np.linalg.inv(ABD2) # Then, calculate abd
Efx2 = 12/(abd2[3][3]*h**3)
Efy2 = 12/(abd2[4][4]*h**3)

# compare the values
print("The found bending stiffness Efs [MPa]: ", Efs)
print("The calculated Efx,1, Efy,1 : ", Efx, ", ", Efy)
print("The calculated EFx,2, Efy,2 : ",Efx2, ", ", Efy2)


# %% A1. layup 1, Fx Max failure load
# Get NM, first determine for layup 1 with Mx
reload(ft)
F = 656
Mx = F*Ls/(4*b_s)# Mx = F L / 4 b
NM1x = [0, 0, 0, Mx, 0, 0]# set NM for this load case
sig_star = cp.CalculateStress(NM1x,C_star_laminate,z) # Use function stresses to get sigma*
sig = cp.RotateMaterial(sig_star,layup1)# rotate stresses to material CS
fails =  ft.MaxStressTest(NM1x, layup1, strength, C_star_laminate, cd['t']) # determine failure ply and fail tests
n_lst = np.array(range(n))
fail_ply_1x = n_lst[fails==True]
fail_ply_1x += 1 
print("The ply that fails is number: ", fail_ply_1x)

# %% A2. Determine failure ply for layup 1 My
# calculate NM1y
F1y = 613
M1y = F1y*Ls/(4*b_s)# Mx = F L / 4 b
NM1y = [0, 0, 0, 0, M1y, 0]# set NM for this load case
sig_star_1y = cp.CalculateStress(NM1y,C_star_laminate,z) # Use function stresses with C*lam1 to get sigma*
sig_1y = cp.RotateMaterial(sig_star_1y,layup1)# rotate stresses to material CS
print("Tested force (F1y): ",str(F1y), "N")
fails_1y =  ft.MaxStressTest(NM1y, layup1, strength, C_star_laminate, cd['t']) # determine failure ply and fail tests
fail_ply_1y = n_lst[fails_1y==True]
fail_ply_1y += 1 
print("The ply that fails is number: ", fail_ply_1y)

# %% A3. Determine failure ply for layup 2 Mx
# calculate NM2x
F2x = 961 # N
M2x = F2x*Ls/(4*b_s)# Mx = F L / 4 b
NM2x = [0, 0, 0, M2x, 0, 0]# set NM for this load case
sig_star_2x = cp.CalculateStress(NM2x,Cstar_lam2,z) # Use function stresses with C*lam1 to get sigma*
sig_2x = cp.RotateMaterial(sig_star_2x,layup2)# rotate stresses to material CS
print("Tested force (F2x): ",str(F2x), "N")
fails_2x =  ft.MaxStressTest(NM2x, layup2, strength, Cstar_lam2, cd['t']) # determine failure ply and fail tests
fail_ply_2x = n_lst[fails_2x==True]
fail_ply_2x += 1 
print("The ply that fails is number: ", fail_ply_2x)

# %% A4. Determine failure ply for layup 2 My
# calculate NM2y
F2y = 439
M2y = F2y*Ls/(4*b_s) # My = F L / 4 b
NM2y = [0, 0, 0, 0, M2y, 0]# set NM for this load case
sig_star_2y = cp.CalculateStress(NM2y,Cstar_lam2,z) # Use function stresses with C*lam1 to get sigma*
sig_2y = cp.RotateMaterial(sig_star_2y,layup2)# rotate stresses to material CS
print("Tested force (F2y): ",str(F2y), "N")
fails_2y =  ft.MaxStressTest(NM2y, layup2, strength, Cstar_lam2, cd['t']) # determine failure ply and fail tests
fail_ply_2y = n_lst[fails_2y==True]
fail_ply_2y += 1 
print("The ply that fails is number: ", fail_ply_2y)

# %% R1. Report on found values layup1, Mx
f = open("FailureTests3PB.txt", "r+")# create new txt file, erase contents
f.truncate(0) # need '0' when using r+

f.write("tested layup : layup 1, moment : Mx")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply_1x)+"\n\n")# Record which ply fails first
f.close()# close file

# %% R2. Report on found values layup1, My
f = open("FailureTests3PB.txt", "a")# append to txt file
f.write("tested layup : layup 1, moment : My")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F1y)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply_1y)+"\n\n")# Record which ply fails first
f.close()# close file

# %% R3. Report on found values layup2, Mx
f = open("FailureTests3PB.txt", "a")# append
f.write("tested layup : layup 2, moment : Mx")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F2x)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply_2x)+"\n\n")# Record which ply fails first
f.close()# close file

# %% R4. Report on found values layup2, My
f = open("FailureTests3PB.txt", "a")# append
f.write("tested layup : layup 2, moment : My")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F2y)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply_2y)+"\n\n")# Record which ply fails first
f.close()# close file

# %% plot stress in longitudinal direction (layup 2, y direction)
reload(cp)
cp.PlotPlyStress(NM2y,Cstar_lam2,z,1)

# %% Calculate thermal stresses
reload(te)
AlphaStarArray = te.alpha_star_laminate(AlphaVec,layup2) # use of 2nd layup
NM_th = te.calculate_NM_thermal(Cstar_lam2, delta, AlphaStarArray, z)
SigTh = te.ThermalStress(NM_th, AlphaStarArray, delta, Cstar_lam2, z) # stresses in ply CS
SigThMat = cp.RotateMaterial(SigTh, layup2) # rotate thermal stresses to material CS
SigThMat = np.round(SigThMat,7)

#PlotThermalStress(SigTh,0,z)
te.PlotThermalStress(SigThMat,0,z)
te.PlotThermalStress(SigThMat,1,z)
# %%
