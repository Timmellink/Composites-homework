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
    if (-s1c<sig1 and sig1<s1t):
        #print("No failure in 1 direction")
        if (-s2c<sig2 and sig2<s2t):
            #print("No failure in 2 direction")
            if sig6<s6:
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

