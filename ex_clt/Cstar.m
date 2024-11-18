function C_star = Cstar(C, theta)
R=  [1 0 0; 0 1 0;0 0 2]; % Reuter matrix
T = transformation(theta);
C_star = ((T\C)*R*T)/R;
end