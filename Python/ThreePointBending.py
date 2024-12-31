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
import TensileTest as tt

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

# %% max stress test 
def MaxStressTest(sig,strength):
    """
    Check per ply whether it fails on MaxStress criterion

    Parameters
    ---------
    sig : matrix
      2xn matrix of stresses in material CS at start and end per layer
    strength : tuple
      A tuple containing s1c,s1t,s2c,s2t and s6
      strength : (S1c,S1t,S2c,S2t,S6)
 
    Returns
    -------
    fail_ply : scalar
      number of last ply in laminate that fails
    Failures : array
      an array of booleans of whether each ply fails
    """
    s1c,s1t,s2c,s2t,s6 = strength
    Failures = []
    for PlyStress in sig:
        sig1 = max(abs(PlyStress[0][0]),abs(PlyStress[1][0])) # get maximum stress in 1-direction of ply
        if -sig1 == min(PlyStress[0][0],PlyStress[1][0]): # value was negative
          sig1 = -sig1
        sig2 = max(abs(PlyStress[0][1]),abs(PlyStress[1][1])) # get max stress in 2 direction of ply
        if -sig2 == min(PlyStress[0][1],PlyStress[1][1]): # value was negative
          sig2 = -sig2
        sig3 = max(abs(PlyStress[0][2]),abs(PlyStress[1][2])) # get stress in 6-direction
        if -sig3 == min(PlyStress[0][2],PlyStress[1][2]): # value was negative
          sig3 = -sig3
        fail = ft.MaxStress(sig1, sig2, sig3, s1c, s1t, s2c, s2t, s6) 
        Failures.append(fail)
    Failures = np.array(Failures) # convert to np array
    n = len(sig) # get number of plies
    n_lst = np.array(range(n)) # get list of ply numbers
    if (len(n_lst[Failures==True])>=1):
      print("The plies that fail are number(s) #", n_lst[Failures==True]) # print which ply fails
      fail_ply = n_lst[Failures==True][-1]
    else:
       print("No ply fails")
       fail_ply = False
    return fail_ply,Failures

# %% A1. layup 1, Fx Max failure load
# Get NM, first determine for layup 1 with Mx
F = 656
Mx = F*Ls/(4*b_s)# Mx = F L / 4 b
NM1x = [0, 0, 0, Mx, 0, 0]# set NM for this load case
sig_star = tt.CalculateStress(NM1x,C_star_laminate,z) # Use function stresses to get sigma*
sig = tt.RotateMaterial(sig_star,layup1)# rotate stresses to material CS
fail_ply,fails =  MaxStressTest(sig,strength) # determine failure ply and fail tests

# %% A2. Determine failure ply for layup 1 My
# calculate NM1y
F1y = 110
M1y = F1y*Ls/(4*b_s)# Mx = F L / 4 b
NM1y = [0, 0, 0, 0, M1y, 0]# set NM for this load case
sig_star_1y = tt.CalculateStress(NM1y,Cstar_lam2,z) # Use function stresses with C*lam1 to get sigma*
sig_1y = tt.RotateMaterial(sig_star_1y,layup1)# rotate stresses to material CS
print("Tested force (F1y): ",str(F1y), "N")
fail_ply_1y,fails_1y =  MaxStress(sig_1y,sd) # determine failure ply and fail tests

# %% A3. Determine failure ply for layup 2 Mx
# calculate NM2x
F2x = 208
M2x = F2x*Ls/(4*b_s)# Mx = F L / 4 b
NM2x = [0, 0, 0, M2x, 0, 0]# set NM for this load case
sig_star_2x = tt.CalculateStress(NM2x,C_star_laminate,z) # Use function stresses with C*lam1 to get sigma*
sig_2x = tt.RotateMaterial(sig_star_2x,layup2)# rotate stresses to material CS
print("Tested force (F2x): ",str(F2x), "N")
fail_ply_2x,fails_2x =  MaxStress(sig_2x,sd) # determine failure ply and fail tests

# %% A4. Determine failure ply for layup 2 My
# calculate NM2y
F2y = 500
M2y = F2y*Ls/(4*b_s) # My = F L / 4 b
NM2y = [0, 0, 0, 0, M2y, 0]# set NM for this load case
sig_star_2y = tt.CalculateStress(NM2y,Cstar_lam2,z) # Use function stresses with C*lam1 to get sigma*
sig_2y = tt.RotateMaterial(sig_star_2y,layup2)# rotate stresses to material CS
print("Tested force (F2y): ",str(F2y), "N")
fail_ply_2y,fails_2y =  MaxStress(sig_2y,sd) # determine failure ply and fail tests


# %% R1. Report on found values layup1, Mx
f = open("FailureTests3PB.txt", "a")# create new txt file
f.write("tested layup : "+str(layup1)+", moment : Mx")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply)+"\n\n")# Record which ply fails first
f.close()# close file

# %% R2. Report on found values layup1, My
f = open("FailureTests3PB.txt", "a")# append to txt file
f.write("tested layup : "+str(layup1)+", moment : My")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F1y)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply_1y)+"\n\n")# Record which ply fails first
f.close()# close file

# %% R3. Report on found values layup2, Mx
f = open("FailureTests3PB.txt", "a")# append
f.write("tested layup : "+str(layup2)+", moment : Mx")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F2x)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply_2x)+"\n\n")# Record which ply fails first
f.close()# close file

# %% R4. Report on found values layup2, My
f = open("FailureTests3PB.txt", "a")# append
f.write("tested layup : "+str(layup2)+", moment : My")# will record tested layup and which moment (Mx or My)
f.write("\nFailure load : "+str(F2y)+"N")# Record failure load
f.write("\nply that failed : "+str(fail_ply_2y)+"\n\n")# Record which ply fails first
f.close()# close file

# %% plot stress in longitudinal direction (y direction)
def PlotStress(PlyStr, z, dir): 
  """
  Plots material stresses of laminate in a certain direction

  Parameters
  ----------
  PlyStr : matrix
    2xn matrix of stresses in ply CS
  z : array
    array of length n+1 of ply edges
  dir : scalar
    direction of desired stresses to plot (1,2 or 3)

  Returns
  -------
  None

  """
  plt.close()
  PlyStrDir = PlyStr[:,:,dir] # stresses in a direction
  PlyStrDirRow = np.reshape(PlyStrDir, (-1))  # get the stresses in one row
  xd = [[z,z]for z in z] # repeat the z coordinates twice
  xdRow = np.reshape(xd,-1) # reshape to one row
  #xdRow =  # reverse order
  plt.plot(PlyStrDirRow,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)
  return 
#PlotStress(sig_star_2y,z,0)
""" 
str2Material = MatStr[:,:,1] # stresses in 2 direction
str2MaterialRow = np.reshape(str2Material, (-1))  # get the stresses in one row
plt.figure()
plt.plot(str2MaterialRow,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)

str3Material = MatStr[:,:,2] # stresses in 3 direction
str3MaterialRow = np.reshape(str3Material, (-1))  # get the stresses in one row
plt.figure()
plt.plot(str3MaterialRow,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)
"""
# %%
