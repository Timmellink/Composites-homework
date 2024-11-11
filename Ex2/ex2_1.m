clc 
clear
%% variables
E1 = 70e9; % GPa
E2 = 70e9;
G12 = 5e9; % GPa
nu12 = 0.25;
sig1_star = 100e6; % MPa
sig2_star = -50e6;
sig6_star = 50e6;
theta = 30; % degrees
%% transform sigma* to sigma
sig_star_vec = [sig1_star; sig2_star; sig6_star]; % vector of sigma_star
T = trnsfrm_matrix(theta);
sig_vec = T*sig_star_vec; % vector of sigma
S = compliance(E1, E2, nu12, G12);
eps_vec = S*sig_vec;

