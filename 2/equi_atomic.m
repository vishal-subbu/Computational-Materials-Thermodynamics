clear all
clc
clf

no_of_comp = 4;
mole_frac = zeros(no_of_comp,1);
mole_frac(1) = 0.25; % mole fraction of Cobalt
mole_frac(2) = 0.25; % mole fraction of Copper
mole_frac(3) = 0.25; % mole fraction of Iron
mole_frac(4) = compute_x4(mole_frac(1),mole_frac(2),mole_frac(3)); % mole fraction of Nickel


del_H_muggianu = 0.0;
del_H_kohler   = 0.0;
del_H_colinet  = 0.0;

% Muggianu Scheme
beta  = 1.0;
for m = 1:no_of_comp-1
    for n = m+1 : no_of_comp
        x_i = ( 1 + mole_frac(m) - mole_frac(n))/2;
        x_j = ( 1 + mole_frac(n) - mole_frac(m))/2;
        f = mole_frac(m)*mole_frac(n)/(x_i*x_j);
        del_H_muggianu = del_H_muggianu + func(m,n,x_i,x_j)*f*beta;
    end
end
% Kohler Scheme
beta  = 1.0;
for m = 1:no_of_comp-1
    for n = m+1 : no_of_comp
        lamda = (mole_frac(m) - mole_frac(n))/(mole_frac(m) + mole_frac(n));
        x_i = (( 1 + mole_frac(m) - mole_frac(n)) + lamda*( 1 - mole_frac(m) - mole_frac(n)))/2;
        x_j = (( 1 + mole_frac(n) - mole_frac(m)) + lamda*( 1 - mole_frac(n) - mole_frac(m)))/2;
        f = mole_frac(m)*mole_frac(n)/(x_i*x_j);
        del_H_kohler= del_H_kohler + func(m,n,x_i,x_j)*f*beta;
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
            del_H_colinet = del_H_colinet + func(m,n,x_i,x_j)*f*beta_colinet(o);
        end
    end
end

% displaying ther output
X = ['Composition '];
disp(X);
X = ['    Co=',num2str(mole_frac(1))];
disp(X);
X = ['    Cu=',num2str(mole_frac(2))];
disp(X);
X = ['    Fe=',num2str(mole_frac(3))];
disp(X);
X = ['    Ni=',num2str(mole_frac(4))];
disp(X);
X = ['Excess Enthalpy in J/mol '];
disp(X);
X = ['    through Muggianu scheme = ',num2str(del_H_muggianu)];
disp(X);
X = ['    through Kohler   scheme = ',num2str(del_H_kohler)];
disp(X);
X = ['    through Colinet  scheme = ',num2str(del_H_colinet)];
disp(X);

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