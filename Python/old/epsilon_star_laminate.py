# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:48:31 2024

@author: Timme
"""

import numpy as np
#import transformation as t
import epsilon_star as es
def epsilon_star_laminate(eps,layup):
    """
    

    Parameters
    ----------
    eps : vector
        Vector of deformations in material CS.
    layup : list
        List of layup (angles of plies in laminate). If symmetric, 
        it needs to be flipped and appended to itself.

    Returns
    -------
    epsilon_star_array : array of rotated deformation vectors (in ply CS).

    """
#    n = len(layup)*2
    theta_reverse = layup[::-1] # reverse layup
    laminate = np.append(layup,theta_reverse)
   # R = np.array([[1, 0, 0 ],[0, 1, 0],[0, 0, 2]]) # Reuter matrix
    #Ri = np.linalg.invert(R) # Inverted Reuter matrix
    epsilon_star_array = [es.epsilon_star(eps,theta) for theta in laminate]
    return epsilon_star_array