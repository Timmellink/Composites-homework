function [S_rotate] = rotate_S(S, T_inv, T, R, R_inv)
S_rotate = R*T_inv*R_inv*S*T;
end