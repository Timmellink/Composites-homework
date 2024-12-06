# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 15:14:57 2024

@author: Timme
"""
import numpy as np
import transformation as t
def epsilon_star(eps,thet):
    """
    

    Parameters
    ----------
    eps : vector
        Vector of deformations in all directions (material CS).
    thet : scalar
        Angle of ply in the laminate.

    Returns
    -------
    eps_star: vector
        Vector of rotated deformations in ply CS.

    """
    R = np.array([[1,0,0],[0,1,0],[0,0,2]])
    Ri = np.linalg.inv(R)
    T = t.transformation(thet)
    Ti = np.linalg.inv(T)
    eps_star = R@Ti@Ri@eps
    return eps_star