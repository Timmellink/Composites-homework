clc
clear
%% variables
E1 = 40e9; % Pa
E2= 10e9;
G12 = 3.5e9;
nu12 = 0.25;
theta = 53.1; % degrees
d = 1; % diameter (m)
R = d/2; % radius4
t = 20e-3; % thickness (mm)
eps_star = 0.001;
%% transform strain in ply direction to strain in material direction
T = trnsfrm_matrix(theta);
%eps_star_vec = [eps1_star; eps2_star; eps6_star];
%eps_vec = T*eps_star_vec;
S  = compliance(E1,E2,nu12,G12);
%sigma_vec = S*eps_vec;
%% calculate pressure
%P = sigma_vec(1)*2*t/d;
reuter = [1 0 0; 0 1 0; 0 0 2];
S_star = reuter*inv(T)*inv(reuter)*S*T;
S_star_12 = S_star(1,2);
S_star_11 = S_star(1,1);
P = t*eps_star/(R*(S_star_11/2+S_star_12))