function [S] = compliance(E_1, E_2, v_12, G_12)
v_21 = ((v_12*E_2)/E_1);
a11 = 1/E_1;
a12 = -v_21/E_2;
a21 = -v_12/E_1;
a22 = 1/E_2;
a13 = 0;
a31 = 0;
a23 = 0;
a32 = 0;
a33 = 1/G_12;
S = [a11 a12 a13; a21 a22 a23; a31 a32 a33];
end