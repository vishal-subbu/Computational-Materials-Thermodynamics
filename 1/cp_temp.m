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
Temp = 0:5:1500;
gibbs= 0:10:1500;
cp   = 0:10:1500;
cp_mo= 0:10:1500;
len = size(Temp,2);
Tc   = 485;
for i = 1:len
    temp = Temp(i);
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
        cp_mo(i) = 2*K_ferro*R*(tau^m + tau^(3*m)/3 + tau^(5*m)/5);
    else
        cp_mo(i)      = 2*K_para*R*(tau^(-n) + tau^(-3*n)/3 + tau^(-5*n)/5);
    end
    
    %% Evaluvation of total thermodynamic properties
    
    cp(i)     = cp_mo(i)      + cp_ch;
    
end
figure(1)
clf
hold on
plot(Temp,cp,'--','linewidth',4.0);
xlabel('Temp(K)');
ylabel('Specific Heat(J mfu^{-1} K^{-1})')
title('Specific heat as a function of temperature')
ax = gca ;
set(ax, 'linewidth',2.0);
axis('square');
hold off
grid on


%% Calculation of G at Tc 
f    = 0.28;
beta = 1.008;
R    = 8.314;
A    = (518/1125) + 11692*(1-f)/(15975*f);
g1   = -((1/10) + (1/315) + (1/1500))/A; 
G_mo = R*Tc*log(1+beta)*g1
g2   = (1/6 + 1/135 + 1/600)*474*(1-f)/(497*f);
g2   = g2 + 79/140;
g2   = g2/A;
G2_mo= R*Tc*log(1+beta)*(1-g2)

for i = 1:len
    tau = Temp(i)/Tc;
    if(tau<=1)
        gibbs(i) = -K_para*R*Tc*(79/140 - 518*tau/1125);
        gibbs(i) = gibbs(i) - K_ferro*R*Tc*(tau^(4)/6 + tau^(10)/135 + tau^(16)/600);
        gibbs(i) = gibbs(i) - K_ferro*R*Tc*(71/120 - 518*tau/675);
    else
        gibbs(i) = -K_para*R*Tc*((tau^(-4)/10) + (tau^(-14)/315) + (tau^(-24)/1500)); 
    end
    if(tau==1)
        gibbs(i)
    end
end
tau = 1;
gibbs_energy = -K_para*R*Tc*((tau^(-4)/10) + (tau^(-14)/315) + (tau^(-24)/1500));
X = ['Gibbs_mo = ',num2str(gibbs_energy),' J mfu^{-1}'];
disp(X);
figure(2)
clf
hold on
plot(Temp,gibbs,'linewidth',4.0);
xlabel('Temp(K)');
ylabel('Gibbs')
title('Gibbs as a function of temperature')
ax = gca ;
set(ax, 'linewidth',2.0);
axis('square');
hold off
grid on


