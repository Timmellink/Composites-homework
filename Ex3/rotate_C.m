function [C_rotate] = rotate_C(C, T_inv, T, R, R_inv)
C_rotate = T_inv*C*R*T*R_inv;
end