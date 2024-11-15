clc
clear
%% params
E1 = 128e9;
E2 = 9e9;
nu12 = 0.35;
G12 = 5e9;
t = 0.15e-3; % thickness layer
layup = [0,45,90,-45]; % symmetric
kx = 0.01; % 1/m
%% calculate ABD
% first, calculate ply edges,
% then, create cell array of rotated stiffness matrices C*
% Finally, calculate ABD, using the first two
n = length(layup)*2; % number of layers
z = ply_edges(t,n); % ply edges array
C_array = Cstar_laminate(layup,n,E1,E2,nu12,G12); % rotated matrices cell array
ABD = ABD_matrix(C_array,z);
%% calculate abd
abd = inv(ABD);
%% calculate N and M, all deformations are zero
% from the calculations in my notebook
%Nx = kx/abd(1,4)/(1-abd(4,4)*abd(1,1)/(abd(1,4)^2));
%Mx = -abd(1,1)/abd(1,4)*Nx;
%% calculate N and M, using submatrices
k0 = [kx; 0; 0];
% extract submatrices from abd matrix
b = abd(1:3,4:6);
a = abd(1:3,1:3);
d = abd(4:6,4:6);
M = ((b\d)-(a\b))\(b\k0);
N = -(a\b)*M;
%% Calculate sigma1* in first layer
% at first layer eps*1 can be calculated
% eps1* = eps0 + z*k
% where z = -h/1+i*t, where i = 1
% Then calculate sig*1 with sig*1 = [C*}1*{eps*1}
h = n*t; % height composite
i_pos = 1; % layer number
z_pos = -h/2 + i_pos*t;
eps_0 = [0; 0; 0];
eps_st1 = eps_0 + z_pos*k0;
sig_st1 = C_array{1}*eps_st1;