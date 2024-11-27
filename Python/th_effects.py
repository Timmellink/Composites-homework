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
    M = np.zeros(3)
    N = np.zeros(3)
    if len(z) -1 == len(alpha_3):
        N = [N+]