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
syms a b c d_1 d_2 d_4 temp
enthalpy_ch = a - c*temp - 3*d_4*(temp^4) - d_2*temp^2 + 2*d_1*temp^(-1);
entropy_ch  = -b -c*(1 + log(temp)) - 4*d_4*temp^3 - 2*d_2*temp + d_1*temp^(-2);
cp_ch       = -c -12*d_4*temp^3 - 2*d_2*temp -2*d_1*temp^(-2);

latex(enthalpy_ch)
latex(entropy_ch)
latex(cp_ch)
