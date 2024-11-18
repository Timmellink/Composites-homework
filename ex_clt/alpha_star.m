function alpha_vec = alpha_star(alph,thet)
R = [1 0 0; 0 1 0; 0 0 2];
T = transformation(thet);
alpha_vec = (R\T)\R*alph;