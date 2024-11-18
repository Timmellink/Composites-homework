import numpy as np
import Cstar as cs
import stiffness as s
# parameters
E1 = 128e9;
E2 = 9e9;
nu12 = 0.35;
G12 = 5e9;
theta = 45;

C = s.stiffness(E1,E2,nu12,G12)
Cst = cs.Cstar(C,theta)
A = np.array([[2, 4], [5, -6]])
B = np.array([[9, -3], [3, 6]])
#C = A + B      # element wise addition
#print(C)

'''
Output:
[[11  1]
 [ 8  0]]
 '''

