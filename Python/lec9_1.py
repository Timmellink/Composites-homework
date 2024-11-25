# %% install packages
import numpy as np
#import matplotlib.pyplot as plt
import ply_edges as pl
import Cstar_laminate as csl
import stiffness as s
import ABD as ab
# %% params
layup = [0,0,90,90]
E1 = 128e9
E2 = 9e9
G12 = 5e9
nu12 = 0.35
t = 0.15e-3
n = len(layup)*2 # number of plies
h = t*n # height of laminate
# %% To calculate ABD, calculate C for ply, then calculate C* array
# Next, calculate z coordinates. Then, calculate ABD by ABD function
C = s.stiffness(E1,E2,nu12,G12)
Cr = csl.Cstar_laminate(layup,C,n)
z = pl.ply_edges(t,n)
ABD = ab.ABD_matrix(Cr,z)
ABD = np.round(ABD,4)
# %% calculate abd, then calculate Ex, Ey, GXY using formulae
abd =  np.linalg.inv(ABD)
Ex = 1/(abd[0][0]*h) # Ex = 1/(a11*h) h is laminate thickness
Ey = 1/(abd[1][1]*h) # Ey = 1/(a22*h)
Gxy = 1/(abd[2][2]*h) # Gxy = 1/(a33*h)
# %% flexular moduli 
Efx = 12/(abd[3][3]*h**3) # Efx = 12/(d11*h^3)
Efy = 12/(abd[4][4]*h**3) # Efy = 12/(d22*h^3)
# %% print values
print("Ex = "+str(Ex)+'\n')
print("Ey = "+str(Ey)+'\n')
print("Gxy = "+str(Gxy)+'\n\n')
print("Efx = "+str(Efx)+'\n')
print("Efy = "+str(Efy)+'\n')