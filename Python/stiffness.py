import numpy as np
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