    
import numpy as np
def thermal_effects(Cst,delta,alpha_r,z):
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