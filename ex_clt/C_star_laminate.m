function array = c_star_laminate(theta_matrix,n,E1,E2,nu12,G12)
% return cell array of c* matrices

C = stiffness(E1,E2,nu12,G12)
array = cell(1,n)
% first calculate C*s for bottom half (from z =-h/2 to 0)
for i = theta_matrix
    C_st = Cstar(C,theta_matrix(i))
