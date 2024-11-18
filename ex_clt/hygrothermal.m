function array = hygrothermal(alpha,theta_matrix,n)
array = cell(1,n); % set up cell array
%% fill array with alpha*s
% first calculate alpha*s for top half (from z =-h/2 to 0)
for i = 1:length(theta_matrix)
    alpha_st = alpha_star(alpha,theta_matrix(i));
    array{i} = alpha_st;
end
% next, calculate alpha*s for bottom half (z: <0,h/2>)
% first, flip theta_matrix
theta_matrix_flip = flip(theta_matrix);
k = length(theta_matrix); % start from i = end of last position
for i = 1:length(theta_matrix_flip)
    alph_st = alpha_star(alpha,theta_matrix_flip(i));
    array{k+i} = alph_st; % add alpha* to cell array
end
end
