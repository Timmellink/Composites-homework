function [T] = transformation(theta)
m = cosd(theta);
n = sind(theta);
a11 = m^2;
a12 = n^2;
a13 = 2*m*n;
a21 = n^2;
a22 = m^2;
a23 = -2*m*n;
a31 = -m*n;
a32 = m*n;
a33 = m^2-n^2;
T = [a11 a12 a13; a21 a22 a23; a31 a32 a33];

end