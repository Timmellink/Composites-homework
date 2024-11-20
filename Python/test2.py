import numpy as np
import matplotlib.pyplot as plt
import ply_edges as pl
import Cstar_laminate as csl
import stiffness as s
import ABD as ab
# %% params
E1 = 128e9
E2 = 9e9
nu12 = 0.35
G12 = 5e9
t = 0.15e-3  # thickness layer
layup = [0, 45, 90, -45]  # symmetric
kx = 0.01  # 1/m
eps_0 = np.array([0, 0, 0])
N_app = np.array([0, 0, 0])
M_app = np.array([0, 200, 0])

# %% calculate ABD
# first, calculate ply edges,
# then, create list of rotated stiffness matrices C*
# Finally, calculate ABD, using the first two
n = len(layup) * 2  # number of layers
z = pl.ply_edges(t, n)  # ply edges array
C = s.stiffness(E1, E2, nu12, G12)
C_array = csl.Cstar_laminate(layup, C, n)  # rotated matrices list
ABD = ab.ABD_matrix(C_array, z)
ABD = np.round(ABD, 4)


# %% calculate abd
abd = np.linalg.inv(ABD)

# %% calculate N and M, all deformations are zero
# calculate N and M, by NM = ABD * eps_k0
k0 = np.array([kx, 0, 0])
eps0 = np.array([0, 0, 0])
eps_k0 = np.concatenate((eps0, k0))
NM = ABD @ eps_k0

# %% Calculate sigma1* in first layer
# at first layer eps*1 can be calculated
# eps1* = eps0 + z*k
# where z = -h/1+i*t, where i = 1
# Then calculate sig*1 with sig*1 = [C*}1*{eps*1}
sig_st = [None] * n  # stress distribution list
h = n * t  # height composite
eps_k = np.linalg.solve(ABD, np.concatenate((N_app, M_app)))
k = eps_k[3:6]  # extracting k vector
for i_pos in range(n):
    z_pos = -h/2 + (i_pos + 1) * t
    eps_st_i = eps_0 + z[i_pos] * k  # epsilon at start layer
    eps_st_i_N = eps_0 + z[i_pos + 1] * k  # epsilon at end layer
    sig_st_i = C_array[i_pos] @ eps_st_i  # stress at start layer
    sig_st_i_N = C_array[i_pos] @ eps_st_i_N  # stress at end layer
    sig_st[i_pos] = (sig_st_i, sig_st_i_N)

# %%plot stress distribution
plt.close('all')
# plot stress in 1* and 2* direction
# first plot stress in 1* direction 
# as a test, scatter plot the 1* direction stress in first layer
# now plot 1* direction for each layer in the same graph
# plot stress from start to end position
x = np.linspace(0, h, n + 1)  # height of laminate subdivided in n+1 elements
sig_st_fl = sig_st[::-1]
for i in range(n):
    plt.plot([sig_st_fl[i][1][0], sig_st_fl[i][0][0]], [abs(x[i]), abs(x[i + 1])], '-bo')
for i in range(n):
    plt.plot([sig_st_fl[i][1][1], sig_st_fl[i][0][1]], [abs(x[i]), abs(x[i + 1])], '-ro')
