# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:19:37 2024

@author: Timme
"""

# module containing all relevant functions for calculation of laminate properties
import numpy as np
import matplotlib.pyplot as plt

def nPlies(layup):
    """
    Returns number of plies for a symmetric layup

    Parameters
    ----------
    layup : array
        Array of layup.

    Returns
    -------
    n : scalar
        Number of plies in laminate.

    """
    return len(layup)*2

def inv(matrix):
    """
    Invert matrix or vector

    Parameters
    ----------
    matrix : matrix
        input matrix.

    Returns
    -------
    matrix
        inverted matrix.

    """
    return np.linalg.inv(matrix)

def ply_edges(z,n):
    """
    ply_edges returns ply top and bottom surface locations
    
    Parameters
    ----------
    z : scalar 
        ply thickness
    n : scalar 
        Number of layers
        
    
    Returns
    -------
    z : array
        Array of length n+1 with locations of ply
    """
    z = np.linspace(-n*z/2,n*z/2,n+1)
    return z



def stiffness(E1,E2,nu12,G12):
    """
    stiffness returns C matrix
    
     Parameters
     ----------
     E1 : scalar 
        stiffness of ply in 1 direction
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
        angle of ply (in degrees)
        
    
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
    T = np.round(T,7)
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
    
     Returns
     -------
     C_star : matrix
        Rotated Stiffness matrix
    """

    #array = [None] * n  # set up array
    #if s:
    #else:  
       # layup = theta_array
    layup = theta_array + theta_array[::-1]
    array  = [Cstar(C,angle) for angle in layup] 
    # fill array with C*s according to angle in symmetric layup
  #  for angle,index in zip(layup,range(n)):
  #      Cstar = Cstar(C, angle)
  #      array[index] = Cstar
    array = np.round(array,7)
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
    ABD = np.round(ABD,9)
    return ABD

# %% calculate stresses in ply CS function
def CalculateStress(NM,Cstar_laminate,z):
    """
    Calculate stresses in plies in ply CS

    Parameters
    ----------
    NM : vector
        Vector of forces and moments per unit width 
    Cstar_laminate : array
        array of C*s
    z : array
        edges of plies

    Returns
    -------
    stresses : array
        (3,n)-dimensional array. Array of stresses in ply CS in 1,2, and 6-direction

    """
    n = len(z)-1 # n is the number of plies
    ABD = ABD_matrix(Cstar_laminate, z)# Calculate ABD matrix
    abd = np.linalg.inv(ABD) # calculate abd    
    epsk0 = abd@NM # Determine eps_k0 based on this load NM, epsk0 = abd@F
    # Get sigma* by determining sig* = C*lam (eps0 + z k0)  
    eps0 = epsk0[0:3] # eps0 = epsk0[0:3]
    k0 = epsk0[3:6] # k0 = epsk0[3:6]
    stresses = []
    # sigma_i_star = C*[i]@(eps0+z[i]*k0)
    for i in range(n):
        sig_k_b = Cstar_laminate[i]@(eps0+z[i]*k0) # start stress
        sig_k_e = Cstar_laminate[i]@(eps0+z[i+1]*k0) # end stress
        stresses.append([sig_k_b,sig_k_e])
    stresses = np.array(stresses)
    return stresses

# %% rotate back to material CS function
def RotateMaterial(sig_star, layup):
    """
    Rotate stresses back to material CS

    Parameters
    ---------
    sig_star : array
        stresses in ply CS in 1,2 and 6 direction
    layup : array
        layup of angles in laminate of one half.

    Returns
    -------
    MaterialStresses : array
        An (2,n) array of the stresses in material CS in 1,2, and 6-direction
    """

    layup_laminate = layup+layup[::-1] # total layup 
    n = len(layup_laminate)
    str_b = sig_star[:,0,:] # stresses at beginning of plies
    str_e = sig_star[:,1,:] # stresses at end of plies
    MaterialStresses = []
    for i in range(n):
        str_m = transformation(layup_laminate[i])@str_b[i] # get start stress in material CS
        str_mEnd = transformation(layup_laminate[i])@str_e[i] # get end stress in material CS
        MaterialStresses.append([str_m,str_mEnd]) # append to list
    MaterialStresses = np.array(MaterialStresses) # convert material stresses to numpy array

    return MaterialStresses

def PlotPlyStress(NM, C_array,z, direc):
    """
    Plot stresses in ply CS in certain direction

    Parameters
    ----------
    NM : array 
        forces on ply per unit width
    C_array : array
        array of rotated Cstar matrices
    z : array
        n+1 array of ply edges
    dir : scalar
        direction to determine ply stress in (0 : x, 1 : y, or 2 : z)

    Returns
    -------
    null
    """
    # create graph space for one graph
    dict1 =  {
        '0' : 'x',
        '1' : 'y',
        '2' : 'z'} # dictionary of directions to plot stress in
    graph, plot1 = plt.subplots(1,1)
    stresses = CalculateStress(NM,C_array,z)
    str1 = stresses[:,:,direc] # stresses in certain direction
    str1Row = np.reshape(str1, (-1))  # get the stresses in one row
    xd = [[z,z]for z in z] # repeat the elements twice
    xdRow = np.reshape(xd,-1) # reshape to one row
    plot1.plot(str1Row,xdRow[1:-1]) #plot from top to bottom (only take begin and end once)
    title = "stress distribution in "+dict1[str(direc)]+" direction"
    plot1.set_title(title)
    plot1.invert_yaxis()
    return 

def PlotMatStr(NM, Cstar_array, z, dir, layup):
    """
    Plot stresses in material CS

    Parameters
    ----------
    NM : vector
        Force vector
    Cstar_array : array
        array of rotated Cstar matrices
    z : array
        ply edges
    dir : scalar
        direction to plot stresses in
    layup : array
        array of angles of layup (one half)
    """
    dict = {
        '0' : '1',
        '1' : '2',
        '2' : '6'
    }
    #plt.close()
    plot = plt.figure()
    stress= CalculateStress(NM,Cstar_array,z)
    MatStr = RotateMaterial(stress,layup)
    str1Material = MatStr[:,:,dir] # stresses in certain direction
    str1MaterialRow = np.reshape(str1Material, (-1))  # get the stresses in one row
    xd = [[z,z]for z in z] # repeat the elements twice
    xdRow = np.reshape(xd,-1) # reshape to one row
    plt.plot(str1MaterialRow,xdRow[1:-1],'r') #plot from top to bottom (only take begin and end once)
    plt.axvline(x = 0, color = 'b', label = 'axvline - full height', linestyle = ':')
    plt.ylim(max(xdRow[1:-1]), min(xdRow[1:-1]))
    plt.title("stress in "+dict[str(dir)]+"-direction")

    return
