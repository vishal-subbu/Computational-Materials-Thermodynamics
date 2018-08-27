syms R Tc temp tau m n f beta K_ferro K_para
cp_mo = 2*K_ferro*R*(tau^m + tau^(3*m)/3 + tau^(5*m)/5);
enthalpy_mo = 2*K_ferro*R*temp*((tau^m)/(m+1) + tau^(3*m)/(3*(3*m+1)) + tau^(5*m)/(5*(5*m+1)));
    enthalpy_mo = -(R*Tc*(71*K_ferro/120 + 79*K_para/140) - enthalpy_mo);
    entropy_mo  = 2*K_ferro*R*((tau^m)/m + tau^(3*m)/(3*(3*m)) + tau^(5*m)/(5*(5*m)));
    entropy_mo  = -(R*log(1+beta) - entropy_mo);
    
    cp_mo       = 2*K_para*R*(tau^(-n) + tau^(-3*n)/3 + tau^(-5*n)/5);