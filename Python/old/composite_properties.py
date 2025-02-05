# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:19:37 2024

@author: Timme
"""

# module containing all relevant functions for calculation of laminate properties
import numpy as np


def ply_edges(h,n):
    """
    ply_edges returns ply top and bottom surface locations
    
    Parameters
    ----------
    h : scalar 
        ply thickness
    n : scalar 
        Number of layers
        
    
    Returns
    -------
    z : array
        Array of length n+1 with locations of ply
    """
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
    """
    stiffness returns C matrix
    
     Parameters
     ----------
     E1 : scalar 
        stifness of ply in 1 direction
     E2 : scalar 
        Stiffness of ply in 2 direction
    nu12 : scalar
        Poisson ratio
    G12 : scalar
        Shear modulus
        
    
     Returns
     -------
     C : matrix
         Stiffness matrix
    """
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
    """
    transformation returns T matrix
    
     Parameters
     ----------
     theta : scalar 
        angle of ply
        
    
     Returns
     -------
     T : matrix
        Transformation matrix
    """
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
    """
    Cstar returns C_star matrix
    
    Parameters
    ----------
    theta : scalar
        Angle of ply
    C : matrix
        Stiffness matrix of ply
        
    
     Returns
     -------
     C_star : matrix
        Rotated Stiffness matrix
    """
    R=  np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]]) # Reuter matrix
    T = transformation(theta)
    Tinv = np.linalg.inv(T)
    Rinv = np.linalg.inv(R)
    C_star = Tinv@C@R@T@Rinv # C* = T^-1*C*R*T*R^-1
    return C_star



def Cstar_laminate(theta_array, C):
    """
    Cstar_laminate returns array of C_star matrices
    
    Parameters
    ----------
    theta_array : array
        Array of layup 
        (expects a symmetrical array, so only one half of the layup)
    C : matrix
        Stiffness matrix of one ply in material CS
    n : scalar
        Number of plies
        
    
     Returns
     -------
     C_star : matrix
        Rotated Stiffness matrix
    """

    #array = [None] * n  # set up array
    layup = theta_array + theta_array[::-1]
    array  = [Cstar(C,angle) for angle in layup] 
    # fill array with C*s according to angle in symmetric layup
  #  for angle,index in zip(layup,range(n)):
  #      Cstar = Cstar(C, angle)
  #      array[index] = Cstar
    array = np.round(array,6)
    return array



def ABD_matrix(Carray, z):
    """
     ABD_matrix returns ABD matrix
    
     Parameters
     ----------
     Carray : list 
         list of stiffness matrices in ply CS
     z   : array 
         array with location of ply edges
    
     Returns
     -------
     ABD : matrix
         ABD matrix
    
    Example
    -------
    >>> ABD_matrix(Cr,Z)
    [[A B], [B D]]
    """
    Carray = np.array(Carray) # convert list to numpy array
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    
    if len(z) - 1 == len(Carray):
        for i in range(len(Carray)):
            A += Carray[i] * (z[i + 1] - z[i])
            B += Carray[i] / 2 * (z[i + 1]**2 - z[i]**2)
            D += Carray[i] / 3 * (z[i + 1]**3 - z[i]**3)
    
    AB = np.concatenate((A,B),axis=1)
    BD = np.concatenate((B,D),axis=1)
    ABD = np.vstack((AB, BD))
    ABD = np.round(ABD,7)
    return ABD
