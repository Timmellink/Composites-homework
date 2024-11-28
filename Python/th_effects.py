# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:25:34 2024

@author: Timme
"""
import numpy as np
def th_effects(Cst,delta,alpha_r,z):
    """
    

    Parameters
    ----------
    Cst : array
        Array of size n of rotated stifness matrices C*s.
    delta : scalar
        Temperature difference.
    alpha_r : array
        array of vectors of rotated thermal coefficients.
    z : array
        Array of size n+1 with location of ply edges.

    Returns
    -------
    NM_th : vector
        N^th and M^th vector (6 by 1).

    """
    #M = np.zeros(3)
    #N = np.zeros(3)
    lst1 = []
    if len(z) -1 == len(alpha_r):
       # for i in alpha_r:
        #    Ni = Cst[i]*alpha_r[i]*(z[i+1]-z[i])
         #   lst1.append(Ni)
        Nx = [Cst[i]*alpha_r[i]*(z[i+1]-z[i]) for i in alpha_r] # calculation of each Nth force of each layer
        Mx = [Cst[i]*alpha_r[i]*(z[i+1]**2-z[i]**2) for i in alpha_r] # calculation of the Mth force of each layer
        Nx = np.array(Nx) # convert to numpy array for faster computation
        Mx = np.array(Mx)
        N = delta*np.sum(Nx) # sum the forces, multiplied by the temperature difference
        M = delta/2*np.sum(Mx)

        #N = [N+]