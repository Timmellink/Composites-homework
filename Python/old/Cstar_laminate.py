import numpy as np
import Cstar as cs

def Cstar_laminate(theta_M, C, n):
    # return list of C* matrices

    # Calculate params 
    # C = stiffness(E1, E2, nu12, G12);
    array = [None] * n  # set up array

    # fill array of C*s
    # first calculate C*s for top half (from z =-h/2 to 0)
    for i in range(len(theta_M)):
        C_st = cs.Cstar(C, theta_M[i])
        array[i] = C_st

    # next, calculate C*s for bottom half (z: <0,h/2>)
    # first, flip theta_matrix
    theta_matrix_flip = theta_M[::-1]
    k = len(theta_M)  # start from i = end of last position
    for i in range(len(theta_matrix_flip)):
        C_st = cs.Cstar(C, theta_matrix_flip[i])
        array[k + i] = C_st  # add C* to list

    return array
