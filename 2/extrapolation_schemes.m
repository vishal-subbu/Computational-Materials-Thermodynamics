
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
clf

no_of_comp = 4;
mole_frac = zeros(no_of_comp,1);
mole_frac(1) = 0.1; % mole fraction of Cobalt
mole_frac(2) = 0.3; % mole fraction of Copper
mole_frac(3) = 0.6; % mole fraction of Iron
mole_frac(4) = compute_x4(mole_frac(1),mole_frac(2),mole_frac(3)); % mole fraction of Nickel

delx = 0:0.001:0.6 ;
len = size(delx);
del_H_muggianu = zeros(len(2),1);
del_H_kohler = zeros(len(2),1);
del_H_colinet = zeros(len(2),1);
for i = 1:len(2)
    mole_frac(3) = 0.6 - delx(i);
    mole_frac(4) = compute_x4(mole_frac(1),mole_frac(2),mole_frac(3));
    X = ['Calculating for Ni =',num2str(mole_frac(4))];
    disp(X);
    % Muggianu Scheme
    beta  = 1.0;
    for m = 1:no_of_comp-1
        for n = m+1 : no_of_comp
            x_i = ( 1 + mole_frac(m) - mole_frac(n))/2;
            x_j = ( 1 + mole_frac(n) - mole_frac(m))/2;
            f = mole_frac(m)*mole_frac(n)/(x_i*x_j);
            del_H_muggianu(i) = del_H_muggianu(i) + func(m,n,x_i,x_j)*f*beta;
        end
    end
    % Kohler Scheme
    beta  = 1.0;
    for m = 1:no_of_comp-1
        for n = m+1 : no_of_comp
            lamda = (mole_frac(m) - mole_frac(n))/(mole_frac(m) + mole_frac(n));
            x_i = (( 1 + mole_frac(m) - mole_frac(n)) + lamda*( 1 - mole_frac(m) - mole_frac(n)))/2;
            lamda = (mole_frac(n) - mole_frac(m))/(mole_frac(m) + mole_frac(n));
            x_j = (( 1 + mole_frac(n) - mole_frac(m)) + lamda*( 1 - mole_frac(n) - mole_frac(m)))/2;
            f = mole_frac(m)*mole_frac(n)/(x_i*x_j);
            del_H_kohler(i) = del_H_kohler(i) + func(m,n,x_i,x_j)*f*beta;
        end
    end
    % Colinet Scheme
    for m = 1:no_of_comp-1
        for n = m+1 : no_of_comp
            lamda_colinet = [1,-1];
            beta_colinet  = [0.5,0.5];
            for o = 1:2
            x_i = (( 1 + mole_frac(m) - mole_frac(n)) + lamda_colinet(o)*( 1 - mole_frac(m) - mole_frac(n)))/2;
            x_j = (( 1 + mole_frac(n) - mole_frac(m)) + lamda_colinet(o)*( 1 - mole_frac(n) - mole_frac(m)))/2;
            f = mole_frac(m)*mole_frac(n)/(x_i*x_j);
            del_H_colinet(i) = del_H_colinet(i) + func(m,n,x_i,x_j)*f*beta_colinet(o);
            end
        end
    end
end
% Plotting 
hold on
plot(delx*100,del_H_muggianu,'r-*','MarkerIndices',1:20:len(2),'linewidth',1.0);
plot(delx*100,del_H_kohler,'b-o','MarkerIndices',1:30:len(2),'linewidth',1.0);
plot(delx*100,del_H_colinet,'g-+','MarkerIndices',1:50:len(2),'linewidth',1.0);
legend('Muggianu Scheme','Kohler Scheme','Colinet Scheme');
xlabel('Composition of Nickel(at. %)');
ylabel('Excess Enthalpy(J/mol)')
title('Excess enthalpy as a function of composition')
ax = gca ;
set(ax, 'linewidth',2.0);
axis('square');
hold off
grid on
hold off

% function to choose the binary
function z = func(m,n,x_i,x_j)
z = 0;
if((m==1)&&(n==2))
    z = CoCu(x_i,x_j);
end
if((m==1)&&(n==3))
    z = CoFe(x_i,x_j);
end
if((m==1)&&(n==4))
    z = CoNi(x_i,x_j);
end
if((m==2)&&(n==3))
    z = CuFe(x_i,x_j);
end
if((m==2)&&(n==4))
    z = CuNi(x_i,x_j);
end
if((m==3)&&(n==4))
    z = FeNi(x_i,x_j);
end
end

% functions to implement different binary systems 
function y = CoCu(x1,x2)
y = x1*x2*(39332 - 1356*(x1-x2) + 7953*(x1-x2)^2 - 1119*(x1-x2)^3);
end

function y = CoFe(x1,x3)
y = x1*x3*(-9312 - 1752*(x1-x3));
end

function y = CoNi(x1,x4)
y = x1*x4*(1331);
end

function y = CuFe(x2,x3)
y = x2*x3*(35626 - 1530*(x2-x3) + 12714*(x2-x3)^2 + 1177*(x2-x3)^3);
end

function y = CuNi(x2,x4)
y = x2*x4*(12049 - 1862*(x2-x4));
end

function y = FeNi(x3,x4)
y = x3*x4*(-18379 + 9228*(x3-x4));
end


% function to calculate the mole fraction of fourth component
function y = compute_x4(x1,x2,x3)
y = 1 - (x1+x2+x3);
end
