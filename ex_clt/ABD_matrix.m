function ABD = ABD_matrix(C_r, z)
    % ABD_matrix returns ABD matrix
    %
    % Arguments:
    %   C_r : cell array of size n with stiffness matrices in ply CS
    %   z   : Array with size n+1 with location of ply edges
    %
    % Returns:
    %   ABD : ABD matrix
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    if length(z) - 1 == length(C_r)
        for i=1:length(C_r)
            A = A + C_r{i}*(z(i+1)-z(i));
            B = B + C_r{i}/2*(z(i+1)^2-z(i)^2);
            D = D + C_r{i}/3*(z(i+1)^3-z(i)^3);
        end
    end
    ABD = [A, B; B, D];
end


