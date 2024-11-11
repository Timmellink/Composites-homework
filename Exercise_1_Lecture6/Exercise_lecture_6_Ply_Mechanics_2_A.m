clear
clc

% A balanced orthotropic ply is made up of 0 and 90 degree fibers woven into
% a fabric and bonded together. The material has a modulus in fiber direction
% E1 = E2 of 70 GPa, a shear modulus G12 of 5 GPa and a Poisson’s ratio
% ν12 of 0.25. Figure 1 shows the fiber orientation and the stress state of the
% material. Calculate the normal strains in the two fiber directions.

%% strategy
% From Sigma* (σ*) to Epsilon (ε), Stress in material cs to strain in ply

%% Starting data
E_1 = 70e9;                         % Stiffness in 1 direction in GPa
E_2 = 70e9;                         % Stiffness in 2 direction in GPa
v_12 = 0.25;                        % Poison ratio in 12 direction
v_21 = v_12 * (E_2/E_1);            % Poison ratio in 21 direction
G_12 = 5e9;                         % Shear modulus in GPa

%% Input stresses
Sigma_Mat_1 = 100e6;                % Stress in Megapascal
Sigma_Mat_2 = -50e6;                % Stress in Megapascal
Sigma_Mat_6 = 50e6;                 % Stress in Megapascal
Sigma_Mat = [Sigma_Mat_1; Sigma_Mat_2; Sigma_Mat_6];

%% Set up the Reuter Matrix and Inverse
R = [1 0 0; 0 1 0; 0 0 2];          % Reuter matrix 
R_inv = inv(R);                     % Inverse of Reuter matrix

theta = 30;                         % Current angle in degrees
m = cosd(theta);     
n = sind(theta);
    
% Transformation matrix for current angle
T = trnsfrm_matrix(m,n);         
T_inv = inv(T);                     % Inverse Transformation matrix

% Calculate the stiffness and compliance matrix
C = stiffness(E_1, E_2, v_12, G_12);   
S = compliance(E_1, E_2, v_12, G_12);

% Calculate the rotated stiffness and compliance matrix
C_rotate = rotate_C(C, T_inv, T, R, R_inv);    
S_rotate = rotate_S(S, T_inv, T, R, R_inv);    

%% Formulate to answer 
% First calculate strain* from Stress* 
% Then convert strain* to strain in ply cs
Epsilon_mat = S_rotate * Sigma_Mat;
Epsilon_ply = R * T * inv(R) * Epsilon_mat;

%% Results
disp('Strain in ply')
disp(Epsilon_ply)