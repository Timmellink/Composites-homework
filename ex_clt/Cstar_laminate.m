function array = Cstar_laminate(theta_M,n,E1,E2,nu12,G12)
% return cell array of c* matrices

%% Calculate params 
C = stiffness(E1,E2,nu12,G12);
array = cell(1,n); % set up array

%% fill array of C*s
% first calculate C*s for top half (from z =-h/2 to 0)
for i = 1:length(theta_M)
    C_st = Cstar(C,theta_M(i));
    array{i}=C_st;
end

% next, calculate C*s for bottom half (z: <0,h/2>)
% first, flip theta_matrix
theta_matrix_flip = flip(theta_M);
k = length(theta_M); % loop counter
for i = 1:length(theta_matrix_flip)
    C_st = Cstar(C,theta_matrix_flip(i));
    array{k+i} = C_st;
end
end
