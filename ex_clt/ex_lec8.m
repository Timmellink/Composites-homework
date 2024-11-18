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
eps_0 = [0; 0; 0];
N_app = [0; 0; 0];
M_app = [0; 200; 0];
%% calculate ABD
% first, calculate ply edges,
% then, create cell array of rotated stiffness matrices C*
% Finally, calculate ABD, using the first two
n = length(layup)*2; % number of layers
z = ply_edges(t,n); % ply edges array
C_array = Cstar_laminate(layup,n,E1,E2,nu12,G12); % rotated matrices cell array
ABD = ABD_matrix(C_array,z);
ABD = round(ABD,4);
%% calculate abd
abd = inv(ABD);
%% calculate N and M, all deformations are zero
% from the calculations in my notebook
%Nx = kx/abd(1,4)/(1-abd(4,4)*abd(1,1)/(abd(1,4)^2));
%Mx = -abd(1,1)/abd(1,4)*Nx;
%% calculate N and M, using submatrices
k0 = [kx; 0; 0];
eps0 = [0; 0; 0];
eps_k0 = [eps0;k0];
NM = ABD*eps_k0;
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
sig_st = cell(1,n); % stress distribution cell array
h = n*t; % height composite
eps_k = ABD\[N_app;M_app];
k = eps_k(4:6); % extracting k vector
for i_pos=1:n
    %i_pos = 1; % layer number
    z_pos = -h/2 + i_pos*t;
    eps_st_i = eps_0 + z(i_pos)*k; % epsilon at start layer
    eps_st_i_N = eps_0 + z(i_pos+1)*k; % epsilon at end layer
    sig_st_i = C_array{i_pos}*eps_st_i; % stress at start layer
    sig_st_i_N = C_array{i_pos}*eps_st_i_N; % stress at end layer
    sig_st{i_pos} = {sig_st_i, sig_st_i_N};
end
%% plot stress distribution
close all
% plot stress in 1* and 2* direction
% first plot stress in 1* direction 
% as a test, scatter plot the 1* direction stress in first layer
% now plot 1* direction for each layer in the same graph
% plot stress from start to end position
hold on
x = linspace(0,h,n+1); % height of laminate subdivided in n+1 elements
sig_st_fl = flip(sig_st);
for i=1:n
    plot([sig_st_fl{i}{2}(1), sig_st_fl{i}{1}(1)],[abs(x(i)),abs(x(i+1))],'-bo');
end
for i=1:n
    plot([sig_st_fl{i}{2}(2), sig_st_fl{i}{1}(2)],[abs(x(i)),abs(x(i+1))],'-ro');
end
hold off
