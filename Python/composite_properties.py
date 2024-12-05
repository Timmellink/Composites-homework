# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:19:37 2024

@author: Timme
"""

# module containing all relevant functions for calculation of laminate properties
import numpy as np


def ply_edges(h,n):
    # ply edges returns ply top and bottom surface locations
    #
    # Arguments:
    #  h: ply thickness
    #  n : number of plies
    #
    # Output:
    #  z: array of length n+1 with locations of ply
    z = np.linspace(-n*h/2,n*h/2,n+1)
    return z



def stiffness(E1,E2,nu12,G12):
    nu21 = ((nu12*E2)/E1)
    a11 = E1/(1-(nu12*nu21))
    a12 = nu21*E1/(1-(nu12*nu21))
    a21 = nu21*E1/(1-nu12*nu21)
    a22 = E2/(1-(nu12*nu21))
    a13, a31, a23, a32 = (0, 0, 0, 0)
    a33 = G12
    C = np.array([[a11, a12, a13],[a21, a22, a23],[a31, a32, a33]])
    return C



def transformation(theta):
    theta = theta/180*np.pi # convert to radians
    m = np.cos(theta)
    n = np.sin(theta)
    a11 = m**2
    a12 = n**2
    a13 = 2*m*n
    a21 = n**2
    a22 = m**2
    a23 = -2*m*n
    a31 = -m*n
    a32 = m*n
    a33 = m**2-n**2
    T = np.array([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
    return T    



def Cstar(C, theta):
    R=  np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]]); # Reuter matrix
    T = transformation(theta);
    Tinv = np.linalg.inv(T)
    Rinv = np.linalg.inv(R)
    C_star = np.dot(np.dot(np.dot(np.dot(Tinv,C),R),T),Rinv);
    return C_star



def Cstar_laminate(theta_M, C, n):
    # return list of C* matrices

    # Calculate params 
    # C = stiffness(E1, E2, nu12, G12);
    array = [None] * n  # set up array

    # fill array of C*s
    # first calculate C*s for top half (from z =-h/2 to 0)
    for i in range(len(theta_M)):
        C_st = Cstar(C, theta_M[i])
        array[i] = C_st

    # next, calculate C*s for bottom half (z: <0,h/2>)
    # first, flip theta_matrix
    theta_matrix_flip = theta_M[::-1]
    k = len(theta_M)  # start from i = end of last position
    for i in range(len(theta_matrix_flip)):
        C_st = Cstar(C, theta_matrix_flip[i])
        array[k + i] = C_st  # add C* to list

    return array
