syms a b c d_1 d_2 d_4 temp
enthalpy_ch = a - c*temp - 3*d_4*(temp^4) - d_2*temp^2 + 2*d_1*temp^(-1);
entropy_ch  = -b -c*(1 + log(temp)) - 4*d_4*temp^3 - 2*d_2*temp + d_1*temp^(-2);
cp_ch       = -c -12*d_4*temp^3 - 2*d_2*temp -2*d_1*temp^(-2);

latex(enthalpy_ch)
latex(entropy_ch)
latex(cp_ch)