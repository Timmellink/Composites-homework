import CompositeProperties as cp
import numpy as np
from importlib import reload

# functions to test on ply failure
def TestStress(val):
    """
    Test whether stress causes failure

    Params
    ------
    val : scalar

    Returns
    -------
    fail : boolean
    """
    if val<1:
        #print("No failure")
        fail = False
    else:
       # print("Failure")
        fail = True
    return fail


def TsaiHill(sig1,sig2,sig6,s1c,s1t,s2c,s2t,s6):
    """
    Test on tsai hill
    Params
    ------
    sig1 : scalar
    sig2 : scalar
    sig6 : scalar
    s1c : scalar
    s1t : scalar
    s2c : scalar
    s2t : scalar
    s6 : scalar

    Returns
    -------
    fail : boolean
    """
    if (sig1<0 and sig2<0): # compressive cases
        val = sig1**2/s1c**2 - sig1*sig2/s1c**2 + sig2**2/s2c**2 + sig6**2/s6**2
        fail = TestStress(val)
    elif (sig1>0 and sig2>0): # tensile cases
        val = sig1**2/s1t**2 - sig1*sig2/s1t**2 + sig2**2/s2t**2 + sig6**2/s6**2
        fail = TestStress(val)
    elif (sig1>0 and sig2<0): # 1 tensile, 2 compressive
        val = sig1**2/s1t**2 - sig1*sig2/s1t**2 + sig2**2/s2c**2 + sig6**2/s6**2
        fail = TestStress(val)
    else: # 1 compressive, 2 tensile
        val = sig1**2/s1c**2 - sig1*sig2/s1c**2 + sig2**2/s2t**2 + sig6**2/s6**2
        fail = TestStress(val)
    return fail
    

def MaxStress(sig1,sig2,sig6,s1c,s1t,s2c,s2t,s6):
    """
    Test on Max Stress criterion
    Params
    ------
    sig1 : scalar
        stress in 1 direction
    sig2 : scalar
        stress in 2 direction
    sig6 : scalar
        shear stress
    s1c : scalar
        compressive strength of fibers
    s1t : scalar
        tensile strength of fibers
    s2c : scalar
        compressive stress in 2 direction
    s2t : scalar
        tensile strength in 2 direction
    s6 : scalar
        shear strength

    Returns
    -------
    fail : boolean
    """
    if (-s1c<sig1 and sig1<s1t):
        #print("No failure in 1 direction")
        if (-s2c<sig2 and sig2<s2t):
            #print("No failure in 2 direction")
            if abs(sig6)<s6:
                #print("No failure in 6 direction")
                fail = False
            else:
                #print("Failure in 6 direction")
                fail = True
        else:
            #print("Failure in 2 direction")
            fail = True
    else:
        #print("Failure in 1 direction")
        fail = True
    return fail

def TsaiTest(NM,layup,strength,Cstar_lam,t):
    """
    Check per ply whether it fails on Tsai-hill test

    Parameters
    ---------
    NM : vector
        NM vector (Newton/unit width)
    layup : array
        layup of angles in laminate of one half.
    strength : tuple
        A 5 dimensional tuple, containing s1c,s1t,s2c,s2t and s6
        strength : (S1c,S1t,S2c,S2t,S6)
    Cstar_lam : array
        rotated stiffness matrix
    t : scalar
        ply thickness

    Returns
    -------
    Failures : array
        an array of booleans of whether each ply fails
    """
    reload(cp)
    n = len(layup)*2
    #h = t*n
    z = cp.ply_edges(t,n)
    s1c,s1t,s2c,s2t,s6 = strength
    stresses = cp.CalculateStress(NM,Cstar_lam,z)
    MaterialStresses = cp.RotateMaterial(stresses,layup)
    Failures = []
    for StressArray in MaterialStresses:
        sig1 = StressArray[0][0] # get stress in 1-direction at start ply
        sig2 = StressArray[0][1] # get stress in 2 direction at start ply
        sig3 = StressArray[0][2] # get stress in 6-direction
        fail = TsaiHill(sig1, sig2, sig3, s1c, s1t, s2c, s2t, s6) 
        Failures.append(fail)
    Failures = np.array(Failures) # convert to np array
    return Failures