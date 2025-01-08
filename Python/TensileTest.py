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



#stresses = CalculateStress(Fx)

# %% plot ply CS stresses
def PlotPlyStress(NM, C_array,z, direc):
    """
    Plot stresses in ply CS in certain direction

    Parameters
    ----------
    NM : array 
        forces on ply per unit width
    C_array : array
        array of rotated Cstar matrices
    z : array
        n+1 array of ply edges
    dir : scalar
        direction to determine ply stress in (0 : x, 1 : y, or 2 : z)

    Returns
    -------
    null
    """
    # create graph space for one graph
    dict1 =  {
        '0' : 'x',
        '1' : 'y',
        '2' : 'z'} # dictionary of directions to plot stress in
    graph, plot1 = plt.subplots(1,1)
    stresses = CalculateStress(NM,C_array,z)
    str1 = stresses[:,:,direc] # stresses in certain direction
    str1Row = np.reshape(str1, (-1))  # get the stresses in one row
    xd = [[z,z]for z in z] # repeat the elements twice
    xdRow = np.reshape(xd,-1) # reshape to one row
    plot1.plot(str1Row,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)
    title = "stress distribution in "+dict1[str(direc)]+" direction"
    plot1.set_title(title)
    plot1.invert_yaxis()
    return 

Fx = 50e3 # Fx = 50kN
Nx = Fx/W # determine NM vector
NM = [Nx, 0, 0, 0, 0, 0]
PlotPlyStress(NM, Cstar_laminate,z,0)




# %% Test whether plies fail

strength = (S1c,S1t,S2c,S2t,S6)
Fx = 50e3 # Fx = 50kN
Nx = Fx/W # determine NM vector
NM = [Nx, 0, 0, 0, 0, 0]
fail = TsaiTest(NM,layup,strength)
n_lst = np.array(range(n))
#Failures,MaterialStresses = RotateMaterial(stresses,layup)
#print("The plies that fail are number(s) #", n_lst[Failures==True])
stress=CalculateStress(NM,Cstar_laminate,z)
MatStr = RotateMaterial(stress,layup)
print("The plies that fail are number(s) #", n_lst[fail==True])

# %% plot number of plies that fail
"""
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
"""
# %% plot material stresses
def PlotMatStr(NM, Cstar_array, z, dir, layup):
    """
    Plot stresses in material CS

    Parameters
    ----------
    NM : vector
        Force vector
    Cstar_array : array
        array of rotated Cstar matrices
    z : array
        ply edges
    dir : scalar
        direction to plot stresses in
    layup : array
        array of angles of layup (one half)
    """
    dict = {
        '0' : '1',
        '1' : '2',
        '2' : '6'
    }
    #plt.close()
    plot = plt.figure()
    stress=CalculateStress(NM,Cstar_array,z)
    MatStr = RotateMaterial(stress,layup)
    str1Material = MatStr[:,:,dir] # stresses in certain direction
    str1MaterialRow = np.reshape(str1Material, (-1))  # get the stresses in one row
    xd = [[z,z]for z in z] # repeat the elements twice
    xdRow = np.reshape(xd,-1) # reshape to one row
    plt.plot(str1MaterialRow,xdRow[1:-1],'r') #plot from top to bottom (only take begin and end once)
    plt.axvline(x = 0, color = 'b', label = 'axvline - full height', linestyle = ':')
    plt.ylim(max(xdRow[1:-1]), min(xdRow[1:-1]))
    plt.title("stress in "+dict[str(dir)]+"-direction")

strength = (S1c,S1t,S2c,S2t,S6)
Fx = 51.9e3 # Fx = 51.9kN
Nx = Fx/W # determine NM vector
NM = [Nx, 0, 0, 0, 0, 0]
PlotMatStr(NM,Cstar_laminate,z,0,layup)
PlotMatStr(NM,Cstar_laminate,z,1,layup)
PlotMatStr(NM,Cstar_laminate,z,2,layup)

