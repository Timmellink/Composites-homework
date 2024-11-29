# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:34:15 2024

@author: Timme
"""
import numpy as np
#import transformation as t
import alpha_star as ast
def alpha_star_laminate(alph,thet):
    """
    

    Parameters
    ----------
    alph : array
        vector of thermal expansion coeffients in material CS.
    thet : array
        array of layup of laminate.

    Returns
    -------
    alpha_array : array.
        array of rotated thermal expansion coefficient vectors in ply CS

    """
    n = len(thet)*2
    thetr = thet[::-1] # inverted theta array
    thet_lam = np.append(thet,thetr)
    alpha_array = [0]*n # python list
    # first calculate alpha coefficients for entire laminate
    for i in range(len(thet_lam)):
        alpha_array[i] = ast.alpha_star(alph,thet_lam[i])
    return alpha_array