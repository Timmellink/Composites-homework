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
        array of thermal expansion coeffients in material CS.
    thet : angle
        Angle of one ply in ply CS.

    Returns
    -------
    alpha_vec : vector.
        Vector of rotated thermal expansion coefficients in ply CS

    """
    R = np.array([[1, 0, 0 ],[0, 1, 0],[0, 0, 2]])
    T = t.transformation(thet)
    alpha_vec = np.linalg.solve(np.linalg.solve(R,T),R)@alph
    return alpha_vec