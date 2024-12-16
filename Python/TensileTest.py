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
"""
for i in abd:
    print("\n")
    for j in i:
        print(f"{j:01} ", end="")
"""
# %% calculate values engineering constants
Ex = 1/(abd[0][0]*h)# Ex = 1/(a11*h) h is laminate thickness
Ey = 1/(abd[1][1]*h)# Ey = 1/(a22*h)
Gxy = 1/(abd[2][2]*h)# Gxy = 1/(a33*h)
nuxy = -(abd[1][0])/(abd[0][0]) # nuxy = -a21/a11
"""
print("Ex : ",str(Ex))
print("Ey : ",str(Ey))
print("Gxy : ",str(Gxy))
print("nu_xy : ",str(nuxy))
"""

# %% calculate stresses in ply CS function
Fx = 50e3 # Fx = 50kN
def CalculateStress(NM,Cstar_laminate,z):
    """
    Calculate stresses in plies in ply CS

    Parameters
    ----------
    NM : vector
        Vector of forces and moments per unit width 
    Cstar_laminate : array
        array of C*s
    z : array
        edges of plies
    abd : matrix
        inverse ABD matrix

    Returns
    -------
    stresses : array
        (3,n)-dimensional array. Array of stresses in ply CS in 1,2, and 6-direction

    """
    n = len(z)-1 # n is the number of plies
    ABD = cp.ABD_matrix(Cstar_laminate, z)# Calculate ABD matrix
    abd = np.linalg.inv(ABD) # calculate abd    
    epsk0 = abd@NM # Determine eps_k0 based on this load NM, epsk0 = abd@F
    # Get sigma* by determining sig* = C*lam (eps0 + z k0)  
    eps0 = epsk0[0:3] # eps0 = epsk0[0:3]
    k0 = epsk0[3:6] # k0 = epsk0[3:6]
    stresses = []
    # sigma_i_star = C*[i]@(eps0+z[i]*k0)
    for i in range(n):
        sig_k_b = Cstar_laminate[i]@(eps0+z[i]*k0) # start stress
        sig_k_e = Cstar_laminate[i]@(eps0+z[i+1]*k0) # end stress
        stresses.append([sig_k_b,sig_k_e])
    stresses = np.array(stresses)
    return stresses

#stresses = CalculateStress(Fx)

# %% plot stresses
"""
str1 = stresses[:,:,0] # stresses in 1 direction
str1Row = np.reshape(str1, (-1))  # get the stresses in one row
xd = [[z,z]for z in z] # repeat the elements twice
xdRow = np.reshape(xd,-1) # reshape to one row
plt.plot(str1Row,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)
"""

# %% rotate back to material CS function
def RotateMaterial(sig_star, layup):
    """
    Rotate stresses back to material CS

    Parameters
    ---------
    sig_star : array
        stresses in ply CS in 1,2 and 6 direction
    layup : array
        layup of angles in laminate of one half.

    Returns
    -------
    MaterialStresses : array
        An (2,n) array of the stresses in material CS in 1,2, and 6-direction
    """

    layup_laminate = layup+layup[::-1] # total layup 
    str_b = sig_star[:,0,:] # stresses at beginning of plies
    str_e = sig_star[:,1,:] # stresses at end of plies
    MaterialStresses = []
    for i in range(n):
        str_m = cp.transformation(layup_laminate[i])@str_b[i] # get start stress in material CS
        str_mEnd = cp.transformation(layup_laminate[i])@str_e[i] # get end stress in material CS
        MaterialStresses.append([str_m,str_mEnd]) # append to list
    MaterialStresses = np.array(MaterialStresses) # convert material stresses to numpy array

    return MaterialStresses

# %% Run Tsai-Hill test function 
def TsaiTest(Fx,layup,strength):
    """
    Check per ply whether it fails on Tsai-hill test

    Parameters
    ---------
    Fx : scalar
        Force applied in 1-direction (Newton)
    layup : array
        layup of angles in laminate of one half.
    strength : tuple
        A 5 dimensional tuple, containing s1c,s1t,s2c,s2t and s6
        strength : (S1c,S1t,S2c,S25,S6)
    Returns
    -------
    Failures : array
        an array of booleans of whether each ply fails
    """
    
    s1c,s1t,s2c,s2t,s6 = strength
    stresses = CalculateStress(Fx,Cstar_laminate,z,abd,W)
    MaterialStresses = RotateMaterial(stresses,layup)
    Failures = []
    for StressArray in MaterialStresses:
        sig1 = StressArray[0][0] # get stress in 1-direction at start ply
        sig2 = StressArray[0][1] # get stress in 2 direction at start ply
        sig3 = StressArray[0][2] # get stress in 6-direction
        fail = ft.TsaiHill(sig1, sig2, sig3, s1c, s1t, s2c, s2t, s6) 
        Failures.append(fail)
    Failures = np.array(Failures) # convert to np array
    return Failures
# %% Test whether plies fail
""" 
strength = (S1c,S1t,S2c,S2t,S6)
fail = RunTest(Fx,layup,strength)
n_lst = np.array(range(n))
#Failures,MaterialStresses = RotateMaterial(stresses,layup)
#print("The plies that fail are number(s) #", n_lst[Failures==True])
str=CalculateStress(Fx,Cstar_laminate,z,abd,W)
MatStr = RotateMaterial(str,layup)
print("The plies that fail are number(s) #", n_lst[fail==True])
"""
# %% plot number of plies that fail
def PlyBreak(leftBound,rightBound,n,strength,layup):
    n_lst = np.array(range(n))
    force_range = np.linspace(leftBound,rightBound,n)
    n_fail_lst = []
    for force in force_range:
        fail = RunTest(force,layup,strength)
        fail_true = n_lst[fail==True]
        n_fail = len(fail_true)
        n_fail_lst.append(n_fail)
    plt.plot(force_range,n_fail_lst)
    return 
# %% plot stresses
""" 
plt.close()
str1Material = MatStr[:,:,0] # stresses in 1 direction
str1MaterialRow = np.reshape(str1Material, (-1))  # get the stresses in one row
xd = [[z,z]for z in z] # repeat the elements twice
xdRow = np.reshape(xd,-1) # reshape to one row
plt.plot(str1MaterialRow,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)

str2Material = MatStr[:,:,1] # stresses in 2 direction
str2MaterialRow = np.reshape(str2Material, (-1))  # get the stresses in one row
plt.figure()
plt.plot(str2MaterialRow,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)

str3Material = MatStr[:,:,2] # stresses in 3 direction
str3MaterialRow = np.reshape(str3Material, (-1))  # get the stresses in one row
plt.figure()
plt.plot(str3MaterialRow,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)
""" 
