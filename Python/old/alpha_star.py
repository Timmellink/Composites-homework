# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:25:34 2024

@author: Timme
"""
import numpy as np
import transformation as t
def alpha_star(alpha_vec,theta):
    """
    

    Parameters
    ----------
    alpha_vec : array
        Array of thermal coefficients in material CS.
    theta : scalar
        Angle of particular ply.

    Returns
    -------
    alpha_star_vec : vector
        Vector of rotated thermal coefficients (ply CS)

    """
    R = np.array([[1, 0, 0 ],[0, 1, 0],[0, 0, 2]])
    Ri = np.linalg.inv(R) # inverted R matrix
    T = t.transformation(theta)
    Ti = np.linalg.inv(T)
    alpha_star_vec =  R@Ti@Ri@alpha_vec # R*T^-1*R^-1*alpha
    return alpha_star_vec
    