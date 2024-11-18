import numpy as np
import stiffness as s
import transformation as t
def Cstar(C, theta):
    R=  np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]]); # Reuter matrix
    T = t.transformation(theta);
    Tinv = np.linalg.inv(T)
    Rinv = np.linalg.inv(R)
    C_star = np.dot(np.dot(np.dot(np.dot(Tinv,C),R),T),Rinv);
    return C_star
