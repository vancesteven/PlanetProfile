function rho_core = coreDensityRuckrieman(p_GPa,T_K,x_0)
% p_GPa = pressure in pascals
% T_K = temperature in Kelvin
% x_0 = mass fraction of sulfur

% conversion to Pascals
p = p_GPa*1e9;

% initial mass fraction of iron sulfide
chi_0 = x_0/0.3647; % Ruckrieman [near (2)]

% first pressure derivative of the isothermal bulk modulus
dK_T0dp = 5.8; % value from Rivoldini 2009 (Table 3 K_T',ref)

% isothermal bulk modulus
K_T0 = (83.4 - 412*x_0 + 592*(x_0).^2)*1e9; % Ruckrieman 2018 (B1)

% density of liquid iron
rho_0_Feliq = 7034.96 - 0.926.*(T_K-1811); % Ruckrieman 2018 (B3)

% density of liquid FeS
rho_0_FeSliq = 4598-0.58.*(T_K-273); % Ruckrieman 2018 (B4)

% reference density
rho_0 = ( (1-chi_0)./rho_0_Feliq + (chi_0)./rho_0_FeSliq )^(-1); % Ruckrieman 2018 (B2)

% pressure scaling factor
s = (1+dK_T0dp.*(p-10^5)./K_T0).^(1./dK_T0dp);

% pressure-dependent density at the core
rho_core = rho_0.*s; % Ruckrieman 2018 (2)