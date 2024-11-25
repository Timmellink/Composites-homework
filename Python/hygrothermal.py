# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:28:15 2024

@author: Timme
"""

import numpy as np
import alpha_star as a
def hygrothermal(alpha, theta_matrix, n):
    """
    

    Parameters
    ----------
    alpha : array
        array of thermal expansion coefficients in 1,2 and 6 direction.
    theta_matrix : array
        Array of angles of plies.
    n : integer
        Number of layers in laminate.

    Returns
    -------
    array : array
        Returns array of arrays of thermal expansion coefficients in ply CS direction.

    """
   # array = np.zeros(n)
    # %% fill array with alpha*s
    # first calculate alpha*s from top half (z=-h/2 to z=0)
    array = [a.alpha_star(alpha,x) for x in theta_matrix ]    
    """
    for i in theta_matrix:
        alpha_st = a.alpha_star(alpha,i)
        array()
    """
    y = theta_matrix[:,-1]
    array2 = [a.alpha_star(alpha,x) for x in y]
    