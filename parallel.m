%% variables
Em = 2; % matrix modulus (GPa)
Ef = 70; % fibre modulus (GPa)
ksi = 1; % ksi variable
ksi2 = 2;
eta = 1;
eta1 = 0.4;
eta2 = 0.6;
eta3 = 0.2;

%% plot modulus in longitudinal direction
%Ec = par(0.5);
close all
vf_m = 0:0.01:1;
Ec_m = par(vf_m,Em,Ef);
plot(vf_m,Ec_m)
%% plot in transverse direction
Ec2_m = tran(vf_m,Em,Ef);
hold on
plot(vf_m,Ec2_m)
hold off
%% plot with halpin-tsai function
figure()
Ec_h = halp(vf_m,Em,Ef,ksi);
Ec_h2 = halp(vf_m,Em,Ef,ksi2);
plot(vf_m,Ec_h)
hold on 
plot(vf_m,Ec_h2)
hold off
ylabel("E-modulus composite (GPa)")
xlabel("fibre volume fraction")
title("Halpin-Tsai")
%% plot with tsai-han function
figure()
%Eta_m = [eta1, eta2, eta3, eta];
Ec_t = han(vf_m,Em,Ef,eta);
plot(vf_m,Ec_t)


%% function definition longitudinal direction
function Ec = par(vf,Em,Ef)
    Ec = Em*(1-vf)+Ef*vf;
end
%% function definition transverse direction
function Ec2 = tran(vf,Em,Ef)
    Ec2 = Em*Ef./(Ef*(1-vf)+vf*Em);
end
%% halpin tsai function
function Ec_h = halp(vf,Em,Ef,ksi)
    eta = ((Ef/Em)-1)./((Ef/Em)+ksi);
    Ec_h = Em*(1+ksi*eta*vf)./(1-eta*vf);
end
%% Tsai han function
function Ec = han(vf,Em,Ef,eta)
    vm = (1-vf);
    Ec = (vf+eta.*vm).*((Ef*Em)./(vf.*Em+eta.*vm.*Ef));
end 
