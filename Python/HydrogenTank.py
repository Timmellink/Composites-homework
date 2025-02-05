# %% Import statements
import numpy as np
import CompositeProperties as cp
from importlib import reload
import FailTest as ft
# %% set properties
# elastic properties
E1 = 130e9
E2 = 10e9
G12 = 4.5e9
nu12 = 0.32
alph1 = 0.16e-6
alph2 = 27.62-6

# strength values
S1t = 1500e6
S1c = 1200e6
S2t = 100e6
S2c = 120e6
S6 = 100e6
strength = [S1c, S1t, S2c, S2t, S6]

epsx0 = 0.0025
layup = [54,-54]
h = 10e-3 # total laminate thickness
n = len(layup)*2
t = h/n # ply thickness
R = 15e-2 
# %% calculate ABD
reload(cp)
C = cp.stiffness(E1,E2,nu12,G12) # calculate stiffness C
z = cp.ply_edges(t,n)# calculate z ply edges
Clam = cp.Cstar_laminate(layup, C)# calculate cstar_lam
ABD = cp.ABD_matrix(Clam, z)# calculate ABD
abd = np.linalg.inv(ABD)# calculate abd
# %% Calculate P
a11 = abd[0][0]
a12 = abd[0][1]
P = epsx0/R*2/(a11+2*a12)# P = epsx0/R*2/(a11+2a12)
print("The pressure is: "+str(P*1e-5)+" bar.")
# %% calculate epsy0
a21 = abd[1][0]
a22 = abd[1][1]
epsy0 = a21*P*R/2+a22*P*R # epsy0 = a21*PR/2+a22*PR
print("The tangential strain is: "+str(epsy0)+" m/m.")
# %% do max stress test
reload(ft)
Pmax = 250e5 # set max pressure
Nx_max = Pmax*R/2 # calculate Nx 
Ny_max = Pmax*R # Calculate Ny 
NM_max = [Nx_max, Ny_max, 0, 0, 0, 0] # set NM
Fails = ft.MaxStressTest(NM_max, layup, strength, Clam, t)
# %% Calculate Clam2
reload(cp)
layup2 = [45,-45]# set layup2
Clam2 = cp.Cstar_laminate(layup2, C)# calculate Clam2
Fails2 = ft.MaxStressTest(NM_max, layup2, strength, Clam2, t)# do fail test with layup2
# %%
