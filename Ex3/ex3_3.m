clc
clear
%% params
E1 = 40e9;
E2 = 10e9;
G12 = 3.5e9;
nu12= 0.25;
s1t = 1500e6;
s1c = 900e6;
s2t = 60e6;
s2c = 80e6;
s6 = 75e6;
D = 1;
r = D/2;
t = 20e-3;
P= 25e5; % max pressure
theta = 85;
%% calculate stress in ply cs, then calculate stress in material cs
sig1_st = P*r/(2*t); % stress in 1 direction
sig2_st = P*r/(t); % stress in 2 direction
sig6_st = 0;
sig_st = [sig1_st; sig2_st; sig6_st]; % stress vector ply
T = trnsfrm_matrix(theta);
% sigma = T * sigma*
sig = T*sig_st;
%% test whether ply has failed
maxStress(sig(1),sig(2),sig(3),s1c,s1t,s2c,s2t,s6);
tsaiHill(sig(1),sig(2),sig(3),s1c,s1t,s2c,s2t,s6);