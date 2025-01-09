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
    if 1 not in Failures:
        print("No failures according to Tsai-hill")
    else:
        print("Failure according to Tsai-hill")
    return Failures

def MaxStressTest(NM, layup, strength, Cstar_lam, t):
    """
    Check per ply whether it fails on maxstress criterion.

    Parameters
    ---------
    NM : vector
        load vector (units : Newton/width for N and Newton for M)
    layup : array
        layup of laminate
    strength : list
        list containing strength values for ply:
        strength : [S1c, S1t, S2c, S2t, S6]
    Cstar_lam  : array
        rotated stiffness matrix
    t : scalar
        Thickness of single ply

    Returns
    -------
    fails : list
        List of boolean values for failures per ply
    """

    n = len(layup)*2
    z = cp.ply_edges(t,n)
    s1c,s1t,s2c,s2t,s6 = strength
    stresses = cp.CalculateStress(NM, Cstar_lam, z)
    sig_mat = cp.RotateMaterial(stresses, layup)
    Fails = []
    for sig_vec in sig_mat:
        # get maximum stress for direction 1
        MaxVal1 = max(abs(sig_vec[0][0]),abs(sig_vec[1][0])) # get maximum stress in 1-direction of ply
        MinVal1 =  min(sig_vec[0][0],sig_vec[1][0]) # get lowest stress value for direction 1
        if MaxVal1 == -MinVal1: # Maximum stress MaxVal1 was negative
            MaxSig1 = -MaxVal1 
        else:
            MaxSig1 = MaxVal1

        # get maximum stress for direction 2
        MaxVal2 = max(abs(sig_vec[0][1]),abs(sig_vec[1][1])) # get max stress in 2 direction of ply
        MinVal2 = min(sig_vec[0][1],sig_vec[1][1])# get lowest stress value for direction 2
        if MaxVal2 == -MinVal2: # value was negative
          MaxSig2 = -MaxVal2
        else:
            MaxSig2 = MaxVal2

        # get maximum stress for direction 6
        MaxVal3 = max(abs(sig_vec[0][2]),abs(sig_vec[1][2])) # get stress in 6-direction
        MinVal2 = min(sig_vec[0][2],sig_vec[1][2])
        if MaxVal3 == -MinVal2: # value was negative
            MaxSig3 = -MaxVal3
        else:
            MaxSig3 = MaxVal3
        fail = MaxStress(MaxSig1, MaxSig2, MaxSig3, s1c, s1t, s2c, s2t, s6) 
        Fails.append(fail)
    if 1 not in Fails:
        print("No ply fails according to Max Stress test.")
    else:
        print("The laminate fails according to Max Stress criterion.")
    return Fails

def MaxStressTestPlyFail(sig,strength):
    """
    Check per ply whether it fails on MaxStress criterion.
    Prints and returns the ply that failed. 

    Parameters
    ---------
    sig : matrix
      2xn matrix of stresses in material CS at start and end per layer
    strength : tuple
      A tuple containing s1c,s1t,s2c,s2t and s6
      strength : (S1c,S1t,S2c,S2t,S6)
 
    Returns
    -------
    fail_ply : scalar
      number of last ply in laminate that fails
    Failures : array
      an array of booleans of whether each ply fails
    """
    s1c,s1t,s2c,s2t,s6 = strength
    Failures = []
    for PlyStress in sig:
        sig1 = max(abs(PlyStress[0][0]),abs(PlyStress[1][0])) # get maximum stress in 1-direction of ply
        if -sig1 == min(PlyStress[0][0],PlyStress[1][0]): # value was negative
          sig1 = -sig1
        sig2 = max(abs(PlyStress[0][1]),abs(PlyStress[1][1])) # get max stress in 2 direction of ply
        if -sig2 == min(PlyStress[0][1],PlyStress[1][1]): # value was negative
          sig2 = -sig2
        sig3 = max(abs(PlyStress[0][2]),abs(PlyStress[1][2])) # get stress in 6-direction
        if -sig3 == min(PlyStress[0][2],PlyStress[1][2]): # value was negative
          sig3 = -sig3
        fail = ft.MaxStress(sig1, sig2, sig3, s1c, s1t, s2c, s2t, s6) 
        Failures.append(fail)
    Failures = np.array(Failures) # convert to np array
    n = len(sig) # get number of plies
    n_lst = np.array(range(n)) # get list of ply numbers
    if (len(n_lst[Failures==True])>=1):
      FailList = n_lst[Failures ==True]
      FailPly = FailList[0] + 1
      print("The ply that fails is number #", FailPly) # print which ply fails
    else:
       print("No ply fails")
       FailPly = False
    return FailPly,Failures
