clear
clc

%% Starting values
E_1 =       40e9;  
E_2 =       10e9;
G_12 =      3.5e9;
v_12 =      0.25;
v_21 =      v_12*(E_2/E_1);

Radius =    0.5;   % meter
h =         0.02;  % meter
theta =     53.1;  % Angle in degrees 

Epsilon_mat = 0.001;  % Strain in material direction 1

%% Standard matrix components
R = [1 0 0; 0 1 0; 0 0 2];          % Reuter matrix 
R_inv = inv(R);                     % Inverse of Reuter matrix

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

%% formulation to answer
% Epsilon_mat = S_rotate(1,1) * ((P*R)/(2*h)) + S_rotate(1,2) * ((P*R)/(h));
% Rewrite this formula to P = .....

P = (Epsilon_mat * h) / (Radius * ((S_rotate(1,1) / 2) + S_rotate(1,2)));
P = P/100000; % Convert from Pa to bar
P = round(P); % Round off to nearest number

disp('The pressure is:')
disp(P)
disp('in bar')