%
% Copyright (c) 2018, Vishal_S
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Title: Computational Materials Thermodynamics
% 
% Developer: Vishal S
% 
% Contact Info: vishalsubbu97@gmail.com
%
clear all
clc

a   = 0.0;
b   = 0.0;
c   = 0.0;
d_4 = 0.0;
d_2 = 0.0;
d_1 = 0.0;
temp = 500;
Tc   = 485;
%% Calculation of the chemical contribution
if((temp>0)&&(temp<=43))
    a   = 11369.937746;
    b   = -5.641259263;
    c   = 0.0;
    d_4 = -8.333*10^(-6);
    d_2 = 0.0;
    d_1 = 0.0;
else
    if ((temp>43)&&(temp<=163))
        a   = 11622.647246;
        b   = -59.537709263;
        c   = 15.74232;
        d_4 = 0.0;
        d_2 = -0.27565;
        d_1 = 0.0;
    else 
        if ((temp>163)&&(temp<=2000))
            a   = -10195.860754;
            b   = 690.949887637;
            c   = -118.47637;
            d_4 = 0.0;
            d_2 = -0.0007;
            d_1 = 590527;
        end
    end
end

enthalpy_ch = a - c*temp - 3*d_4*(temp^4) - d_2*temp^2 + 2*d_1*temp^(-1);
entropy_ch  = -b -c*(1 + log(temp)) - 4*d_4*temp^3 - 2*d_2*temp + d_1*temp^(-2);
cp_ch       = -c -12*d_4*temp^3 - 2*d_2*temp -2*d_1*temp^(-2);
%% Calculation of magnetic contribution

f    = 0.28;
beta = 1.008;
R    = 8.314;
K_para   = -(79875/74)*f*log(1+beta)/(293*f - 790);
K_ferro  = (474/497)*(1-f)*(K_para/f);
m        = 3;
n        = 5;
tau      = temp/Tc;
if (temp<=Tc)
    cp_mo = 2*K_ferro*R*(tau^m + tau^(3*m)/3 + tau^(5*m)/5);
    enthalpy_mo = 2*K_ferro*R*temp*((tau^m)/(m+1) + tau^(3*m)/(3*(3*m+1)) + tau^(5*m)/(5*(5*m+1)));
    enthalpy_mo = -(R*Tc*(71*K_ferro/120 + 79*K_para/140) - enthalpy_mo);
    entropy_mo  = 2*K_ferro*R*((tau^m)/m + tau^(3*m)/(3*(3*m)) + tau^(5*m)/(5*(5*m)));
    entropy_mo  = -(R*log(1+beta) - entropy_mo);
else
    cp_mo       = 2*K_para*R*(tau^(-n) + tau^(-3*n)/3 + tau^(-5*n)/5);
    enthalpy_mo = 2*K_para*R*((temp*tau^(-n) - Tc)/(-n+1) + (temp*tau^(-3*n)-Tc)/(3*(-3*n+1)) + (temp*tau^(-5*n)-Tc)/(5*(-5*n+1)));
    enthalpy_mo = enthalpy_mo + 2*K_ferro*R*Tc*(1/(m+1) + 1/(3*(3*m+1)) + 1/(5*(5*m+1)));
    enthalpy_mo = -(R*Tc*(71*K_ferro/120 + 79*K_para/140) - enthalpy_mo);
    entropy_mo  = 2*K_para*R*((tau^(-n)-1)/(-n) + (tau^(-3*n) -1)/(3*(-3*n)) + (tau^(-5*n)-1)/(5*(-5*n)));
    entropy_mo  = entropy_mo + 2*K_ferro*R*(1/m + 1/(3*(3*m)) + 1/(5*(5*m)));
    entropy_mo  = -(R*log(1+beta) - entropy_mo);
end

%% Evaluvation of total thermodynamic properties

 enthalpy = enthalpy_mo + enthalpy_ch;
% enthalpy_ch
% enthalpy_mo
 entropy  = entropy_mo  + entropy_ch;
% entropy_mo
% entropy_ch
 cp       = cp_mo       + cp_ch;
% cp_mo
% cp_ch

X = ['temp = ',num2str(temp)];
disp(X);
X = ['Enthalpy = ',num2str(enthalpy),' J mfu^{-1}'];
disp(X);
X = ['Entropy  = ',num2str(entropy),' J K^{-1}'];
disp(X);
X = ['Cp       = ',num2str(cp),' J mfu^{-1} K^{-1}'];
disp(X);
