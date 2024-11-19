import numpy as np

def ABD_matrix(C_r, z):
    # ABD_matrix returns ABD matrix
    #
    # Arguments:
    #   C_r : list of stiffness matrices in ply CS
    #   z   : array with location of ply edges
    #
    # Returns:
    #   ABD : ABD matrix
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    
    if len(z) - 1 == len(C_r):
        for i in range(len(C_r)):
            A += C_r[i] * (z[i + 1] - z[i])
            B += C_r[i] / 2 * (z[i + 1]**2 - z[i]**2)
            D += C_r[i] / 3 * (z[i + 1]**3 - z[i]**3)
    
    ABD = np.array([[A, B], [B, D]])
    return ABD
