function T = conductiveMantleTemperature(r,rb,a,kr,rho,T0,Q,qb)
% from TurcotteSchubert1982, as described in Cammarano et al. 2006
% kr = 4; % W/m2
% rho = 3300; %kg/m3
% a = 1452e3; %top of mantle, m
m = rho*4/3*pi*(a^3-rb^3); % mass of the mantle in kg
H = Q/m; % internal tidal heating
T = T0+rho*H/6/kr*(a.^2-r.^2)+(rho.*H.*rb.^3./3./kr-qb.*rb.^2./kr)*(1./a-1./r); 
% unsure about the sign at the end. switching 1/a and 1/r signs provides agreement with Cammarona 2006 cold case.



