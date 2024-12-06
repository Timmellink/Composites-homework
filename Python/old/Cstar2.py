import numpy as np
import transformation as t

def Cstar(C, theta):
    R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])  # Reuter matrix
    T = t.transformation(theta)
    C_star = np.linalg.solve(T, C) @ R @ T @ np.linalg.inv(R)
    return C_star
