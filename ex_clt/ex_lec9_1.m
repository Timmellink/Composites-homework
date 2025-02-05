clc
clear
%% params
layup = [0 0 90 90];
E1 = 128e9;
E2 = 9e9;
G12 = 5e9;
nu12 = 0.35;
t = 0.15e-3;
n = length(layup)*2; % number of plies
h=t*n; % laminate thickness
%% To calculate ABD, calculate C for ply, then calculate C* array
% Next, calculate the z coordinates. Then, calculate ABD by
% using ABD function
C = stiffness(E1,E2,nu12,G12);
C_r = Cstar_laminate(layup,C,n);
z = ply_edges(t,n);
ABD = ABD_matrix(C_r,z);
ABD = round(ABD,4);
%% calculate abd, then calculate Ex,Ey,Gxy using formulae
abd = inv(ABD);
Ex = 1/(abd(1,1)*h); % Ex = 1/(a11*h) h is the laminate thickness
Ey = 1/(abd(2,2)*h); % Ey = 1/(a22*h)
Gxy = 1/(abd(3,3)*h); % Gxy = 1/(a33*h)
%% flexular moduli
Efx = 12/(abd(4,4)*h^3); % Ef,x = 12/(d11*h^3)
Efy = 12/(abd(5,5)*h^3); % Ef,y =  12/(d22*h^3)
