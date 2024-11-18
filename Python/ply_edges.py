import numpy as np
def ply_edges(h,n):
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

