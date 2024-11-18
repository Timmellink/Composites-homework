function NM_th = Th_effects(Cst,delta,alpha_r,z)
    % Th_effects returns N^Th and M^Th
    %
    % Arguments:
    %   alpha_r : cell array of size n with alpha vectors in ply CS
    %   z   : Array with size n+1 with location of ply edges
    %
    % Returns:
    %   NM_th : N^Th and M^th vector (6 by 1)
    N = zeros(3,1);
    M = zeros(3,1);
    if length(z) - 1 == length(alpha_r)
        for i=1:length(alpha_r)
            N = N + delta*Cst{i}*alpha_r{i}*(z(i+1)-z(i));
            M = M + delta*Cst{i}*alpha_r{i}/2*(z(i+1)^2-z(i)^2);
        end
    end
    NM_th = [N; M];
end


