# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:03:00 2025

@author: Timme
"""

# %% import statements
import numpy as np
import CompositeProperties as cp
import FailTest as ft

# %% set properties
layup = [54, -54]
h = 4e-3
d = 50e-2

# material properties
E1 = 120e9
E2 = 9.5e9
G12 = 3.5e9
nu12 = 0.32
alph1 = 1e-6
alph2 = 30e-6

# %% calculate abd
C = cp.stiffness(E1, E2, nu12, G12)
n = len(layup)*2
t = h/n
z = cp.ply_edges(t, n)
Cstar_lam = cp.Cstar_laminate(layup, C)
ABD = cp.ABD_matrix(Cstar_lam, z)
abd = np.linalg.inv(ABD)

# %% determine strain values
P = 50e5
R = d/2
a11 = abd[0][0]
a12 = abd[0][1]
#a21 = abd[1][0]
#a22 = abd[1][1]

sig_x = P*R/(2*h)
sig_theta = P*R/h
Nx = sig_x*h
Ny = sig_theta*h
NM = [Nx, Ny, 0, 0, 0, 0]
eps_k = abd@NM
epsx = eps_k[0]
epsTheta = eps_k[1]
print(f"The value for epsx = {epsx} and for epsTheta = {epsTheta}")

# %% calculate max pressure
A11 = ABD[0][0]
A21 = ABD[1][0]
eps1 = 0.01
#Pmax = 2/(3*R)*eps1*(A11+A21)

Re = np.array([[1,0,0],[0,1,0],[0,0,2]]) 
Rinv = np.linalg.inv(Re)
T = cp.transformation(-45)
K = Re@T@Rinv
K11 = K[0][0]
K12 = K[0][1]
P_max = eps1/(R*(a11*K11/2+a12*K12))

print(f"The max. allowed value for P is : {P_max*1e-5} bar")

sig_x_max = P_max*R/(2*h)
sig_theta_max = P_max*R/h
Nxm = sig_x_max*h
Nym = sig_theta_max*h
N = [Nxm, Nym, 0]
eps = K@N