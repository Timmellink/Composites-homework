clear
clc

%% Define the starting variables
E_1 = 128;                          % Stiffness in 1 direction in GPa
E_2 = 9;                            % Stiffness in 2 direction in GPa
v_12 = 0.35;                        % Poison ratio in 12 direction
v_21 = v_12;                        % Poison ratio in 21 direction
G_12 = 5;                           % Shear modulus in GPa
theta_range = 0:90;                 % Angle range from 0 to 180 degrees

%% Preallocate arrays for rotated values
E_1_rotate = zeros(size(theta_range));
E_2_rotate = zeros(size(theta_range));
G_12_rotate = zeros(size(theta_range));
v_12_rotate = zeros(size(theta_range));

%% Set up the Reuter Matrix and Inverse
R = [1 0 0; 0 1 0; 0 0 2];          % Reuter matrix 
R_inv = inv(R);                     % Inverse of Reuter matrix

%% Calculate the stiffness and compliance matrix at each angle
for i = 1:length(theta_range)
    theta = theta_range(i);          % Current angle in degrees
    m = cosd(theta);     
    n = sind(theta);
    
    % Transformation matrix for current angle
    T = trnsfrm_matrix(m,n);         
    T_inv = inv(T);                 % Inverse Transformation matrix

    % Calculate the stiffness and compliance matrix
    C = stiffness(E_1, E_2, v_12, G_12);   
    S = compliance(E_1, E_2, v_12, G_12);

    % Calculate the rotated stiffness and compliance matrix
    C_rotate = rotate_C(C, T_inv, T, R, R_inv);    
    S_rotate = rotate_S(S, T_inv, T, R, R_inv);    

    % Calculate rotated material properties
    E_1_rotate(i) = 1/S_rotate(1,1);  % Rotated E_1
    E_2_rotate(i) = 1/S_rotate(2,2);  % Rotated E_2
    G_12_rotate(i) = 1/S_rotate(3,3); % Rotated G_12
    v_12_rotate(i) = -S_rotate(1,2)/S_rotate(1,1); % Rotated G_12
end

%% Plot the results
figure;

% Plot E1* and E2* over an angle from 0 to 90 degrees
hold on
plot(theta_range, E_1_rotate, 'r-', 'LineWidth', 1);
plot(theta_range, E_2_rotate, 'b--', 'LineWidth', 1);
hold off

xlabel('Angle (degrees)');
ylabel('Rotated Material Properties (GPa)');
legend('E_1 (Rotated)', 'E_2 (Rotated)', 'Location', 'Best');
title('Rotated Material Properties vs. Angle (0° to 90°)');
grid on;

% Plot G12* over an angle from 0 to 90 degrees
figure; 
plot(theta_range, G_12_rotate, 'g-', 'LineWidth', 1);

xlabel('Angle (degrees)');
ylabel('Rotated Material Properties (GPa)');
legend('G_{12} (Rotated)', 'Location', 'Best');
title('Rotated Material Properties vs. Angle (0° to 90°)');
grid on;

% Plot v12 over an angle from 0 to 90 degrees
figure;
plot(theta_range, v_12_rotate, 'g-', 'LineWidth', 1);

xlabel('Angle (degrees)');
ylabel('Rotated Material Properties (GPa)');
legend('v_{12} (Rotated)', 'Location', 'Best');
title('Rotated Material Properties vs. Angle (0° to 90°)');
grid on;
