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


a1         = 0.5;
a2         = 0.5;
comp_of_al = 0.4;
Temp       = 1200;
R          = 8.314;

syms y1 % y1 = site fraction of Al in first sub lattice
syms y2 % y2 = site fraction of Al in second sub lattice
syms tot_gibbs % function to store the final gibbs lattice

y2 = (comp_of_al*(a1+a2) - a1*y1)/(a2); % Conservation of composition

%% Calculater the Gibbs_energy 
Gibbs_ref =  y1*y2*G_Al_Al(Temp)               + ...
             (y1 + y2 -2*y1*y2)*G_Ni_Al(Temp)  + ...
             (1-y1)*(1-y2)*G_Ni_Ni(Temp); 
% Calculate entropy contribution                
G_config  = R*Temp*(a1*(y1*log(y1)+ (1-y1)*log(1-y1)) + ...
                                  a2*(y2*log(y2)+ (1-y2)*log(1-y2)));

% Calculate excess contribution
G_e       = y2*(1-y2)*y1*(L_0_Al__Al_Ni(Temp)           + ...
                                         L_1_Al__Al_Ni(Temp)*(2*y2-1));
G_e       = G_e + y2*(1-y2)*(1-y1)*(L_0_Ni__Al_Ni(Temp) + ...
                                         L_1_Ni__Al_Ni(Temp)*(2*y2-1));
G_e       = G_e + y1*(1-y1)*y2*(L_0_Al_Ni__Al(Temp)     + ...
                                         L_1_Al_Ni__Al(Temp)*(2*y1-1));
G_e       = G_e + y1*(1-y1)*(1-y2)*(L_0_Al_Ni__Ni(Temp) + ...
                                         L_1_Al_Ni__Ni(Temp)*(2*y1-1));


Gibbs_tot = Gibbs_ref + G_config + G_e;
% Find y1 that gives the minimum Gibbs
tot_gibbs = symfun((diff(Gibbs_tot,y1)),y1);
gibbs_double_diff = diff(tot_gibbs,y1);
tot_gibbs = symfun(simplify(tot_gibbs),y1);
total_gibbs = matlabFunction(tot_gibbs);
options = optimoptions('fsolve','FunctionTolerance',1E-12,'StepTolerance',1E-12);
S = fsolve(total_gibbs,0.000005,options);
gibbs_double_diff = diff(tot_gibbs,y1);
%% Print the Output
y1 = S;
if(subs(gibbs_double_diff)) > 0
    disp ('You have reached a minima');
else
    disp ('You have reached a maxima, change the initial value to reach a minima');
end
y2 = (comp_of_al*(a1+a2) - a1*y1)/(a2);
X = ['Site fraction of Al in 1st sublattice = ',num2str(y1)];
disp(X);
X = ['Site fraction of Al in 2nd sublattice = ',num2str(y2)];
disp(X);

%% Function definitions for different parameters

function y = GHSERAL(T)
y = -11278.4 + 188.684*T - 31.7482*T*log(T) - 1.231e+028*(T^( - 9));
end

function y = GHSERNI(T)
y = -5179.16+117.854*T - 22.096*T*log(T) -0.0048407*T^2;
end

function y = G_Al_Al(T)
y=10083-4.813*T+GHSERAL(T);
end

function y = G_Ni_Ni(T)
y=8715.08-3.556*T + GHSERNI(T);
end

function y =G_Ni_Al(T)
y=-56500-10.7*T +1.4975*T*log(T)+0.5*GHSERAL(T)+0.5*GHSERNI(T);
end



function y = L_0_Al_Ni__Ni(~)
y = -22050;
end
function y = L_1_Al_Ni__Ni(~)
y = 1115;
end
function y = L_0_Al_Ni__Al(T)
y = -14225-5.625*T;
end
function y = L_1_Al_Ni__Al(~)
y = 0.0;
end
function y = L_0_Ni__Al_Ni(~)
y = -22050;
end
function y = L_1_Ni__Al_Ni(~)
y = 1115;
end
function y = L_0_Al__Al_Ni(T)
y = -14225-5.625*T;
end
function y = L_1_Al__Al_Ni(~)
y = 0.0;
end

