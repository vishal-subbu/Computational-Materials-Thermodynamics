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
fun = @gibbs;
y1 = fsolve(fun,0.5)
y2 = (comp_of_al*(a1+a2) - a1*y1)/(a2)
function g = gibbs(y1)
g = (78529695623869481*y1)/1374389534720 - (4157*log(y1 + 1/5))/1000 + (4157*log(1 - y1))/1000 + (4157*log(4/5 - y1))/1000 - (4157*log(y1))/1000 + 21408*y1^2 - 17840*y1^3 - 86375810599679017/3435973836800;
end
