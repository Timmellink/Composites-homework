    
# all relevant functions to calculate termal forces


import numpy as np
from CompositeProperties import transformation as t
import CompositeProperties as cp


def alpha_star(alpha_vec,theta):
    """
    Calculate alpha star array for single ply

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
    T = t(theta)
    Ti = np.linalg.inv(T)
    alpha_star_vec =  R@Ti@Ri@alpha_vec # R*T^-1*R^-1*alpha
    return alpha_star_vec



def alpha_star_laminate(alph,thet):
    """
    Calculate alpha star array for laminate

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
        alpha_array[i] = alpha_star(alph,thet_lam[i])
    return alpha_array



def calculate_NM_thermal(Cst, delta, alpha_r, z):
    """
    Calculate Thermal forces NM^Th

    Parameters
    ---------
    Cst : array
        Rotated stiffness matrix for laminate
    delta : scalar
        Temperature difference 
    alpha_r : array
        alpha star array for laminate
    z : array
        Ply edges

    Returns
    -------
    NM : array
        Force vector for thermal forces
    """
    lst1 = []
    lst2 = []
    if len(z) -1 == len(alpha_r):
       # for i in alpha_r:
        #    Ni = Cst[i]*alpha_r[i]*(z[i+1]-z[i])
         #   lst1.append(Ni)
        for i in range(len(alpha_r)):
            Ni = Cst[i]@alpha_r[i]*(z[i+1]-z[i])
            lst1.append(Ni)
            Mi = Cst[i]@alpha_r[i]*(z[i+1]**2-z[i]**2)
            lst2.append(Mi)
        #Nx = [Cst[i]@alpha_r[i]*(z[i+1]-z[i]) for i in alpha_r] # calculation of each Nth force of each layer
        #Mx = [Cst[i]@alpha_r[i]*(z[i+1]**2-z[i]**2) for i in alpha_r] # calculation of the Mth force of each layer
        #Nx = np.array(Nx) # convert to numpy array for faster computation
        #Mx = np.array(Mx)
        Nx = np.array(lst1)
        Mx = np.array(lst2)
        N = delta*np.sum(Nx,axis=0) # sum the forces, multiplied by the temperature difference
        M = delta/2*np.sum(Mx,axis=0)
        NM = np.append(N,M)
        return NM
        #N = [N+]

def ThermalStress(NMth, alphaR, delta, Cstar, z):
    """
    Calculate thermal stresses (only taking into account thermal forces)

    Parameters
    ----------
    NMth : array
        Array of thermal forces
    alphaR : array
        Rotated thermal coefficients for laminate
    delta : scalar
        temperature difference
    Cstar : array
        rotated stiffness matrix laminate
    z : array
        ply edges

    Returns
    -------
    ThStress : array 
        matrix of thermal stresses for laminate
    """
    ABD = cp.ABD_matrix(Cstar, z)
    abd = np.linalg.inv(ABD)
    EpsK0 = NMth@abd
    Eps0 = EpsK0[0:3] 
    K0 = EpsK0[3:6]
    n = len(alphaR)
    ThStress = []
    for i in range(n):
        ThStress_b = Cstar[i]@(Eps0+z[i]*K0-alphaR[i]*delta)
        ThStress_e = Cstar[i]@(Eps0+z[i+1]*K0-alphaR[i]*delta)
        ThStress.append([ThStress_b, ThStress_e])
    ThStress = np.array(ThStress)
    return ThStress