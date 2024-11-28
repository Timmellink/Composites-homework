# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:34:15 2024

@author: Timme
"""
import numpy as np
import transformation as t
def alpha_star(alph,thet):
    """
    

    Parameters
    ----------
    alph : array
        vector of thermal expansion coeffients in material CS.
    thet : array
        array of layup of laminate.

    Returns
    -------
    alpha_vec : array.
        array of rotated thermal expansion coefficient vectors in ply CS

    """
    n = len(thet)*2
    thetr = thet[::-1] # inverted theta array
    thet_lam = np.append(thet,thetr)
    R = np.array([[1, 0, 0 ],[0, 1, 0],[0, 0, 2]])
    Ri = np.linalg.inv(R) # inverted R matrix
    #alpha_array = [np.linalg.solve(np.linalg.solve(R,t.transformation(thet[i])),R)@alph for i in thet]
    alpha_array = [0]*n # python list
    # first calculate alpha coefficients for entire laminate
    for i in range(len(thet_lam)):
        T = t.transformation(thet_lam[i])
        Ti = np.linalg.inv(T) # rotated transformation matrix
        # rotate alpha vector in material cs to ply CS
        # alpha* = R^-1*T^-1*R*alpha_vec
        alpha_vec = Ri@Ti@R@alph
        alpha_array[i] = alpha_vec
    return alpha_array