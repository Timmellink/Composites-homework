clc
clear
%% params
E1 = 128e9;
E2 = 9e9;
G12 = 5e9;
nu12 = 0.35;
t = 0.15e-3;
alp1 = 2e-6;
alp2 = 50e-6;
layup = [0 45 90 -45];
delta_T = 100;
alpha = [alp1;alp2;0];
n = length(layup)*2;
%% calculate alpha star cell array
% first, calculate C
% then, calculate C star laminate array
% next, calculate alpha cell array, and use delta T, alpha and C star
% array, and z to calculate NM_th
C= stiffness(E1,E2,nu12,G12);
Cstar_lam = Cstar_laminate(layup,C,n);
z = ply_edges(t,n); 
alpha_array = hygrothermal(alpha,layup,n);
NM = Th_effects(Cstar_lam,delta_T,alpha_array,z);
NM = round(NM,7);