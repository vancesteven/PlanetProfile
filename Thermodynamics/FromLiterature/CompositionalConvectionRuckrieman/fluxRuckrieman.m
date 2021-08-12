function [F_l,F_c] = fluxRuckrieman( inputVals )
%{
    Calculates the heat flux out of mantle in core in accordance with
    appendix B2 of Ruckrieman 2018.

    Input:
        Struct that potentially contains fields for independent variables
            T_m (Temperature at upper mantle) [K]
            T_c (Temperature at the CMB) [K]
            T_int (Temperature at outer radius of mantle) [K]
            k_m (Mantle thermal conductivity) [W / m*K]
            E_m (Activation energy of silicate mantle) [J / mol]
            R_m (Outer radius of the mantle) [m]
            R_c (Radius of the core) [m]
            D_l (Thickness of stagnant lid) [m]
            rho_m (Silicate mantle density) [kg / m^3]
            rho_c (Core density) [kg / m^3]
            alpha_m (Mantle thermal expansivity) [1 / K]
            c_pm (Mantle heat capacity) [J / K*kg]
            kappa_m (Mantle thermal diffusivity) [m^2 / s]
            nu (Mantle viscocity) [Pa * s]
%}


% Temperature at upper mantle [K]
if isfield( inputVals , 'T_m' )
    T_m = inputVals.T_m;
else
    T_m = 1750;
end
% Temperature at CMB [K]
if isfield( inputVals , 'T_c' )
    T_c = inputVals.T_c;
else
    T_c = 2000;
end
% Temperature at outer radius of mantle [K]
if isfield( inputVals , 'T_int' )
    T_int = inputVals.T_int;
else
    T_int = 250;
end
% Mantle thermal conductivity [W / m*K]
if isfield( inputVals , 'k_m' )
    k_m = inputVals.k_m;
else
    k_m = 4; % cites Hauck et al. 2006
end
% Activation energy of silicate mantle [J / mol]
if isfield( inputVals , 'E_m' )
    E_m = inputVals.E_m;
else
    E_m = 240*10^3; % cites Karato and Wu 1993
end
% Outer radius of the mantle [m]
if isfield( inputVals , 'R_m' )
    R_m = inputVals.R_m;
else
    R_m = 1770*10^3;
end
% Radius of the core [m]
if isfield( inputVals , 'R_c' )
    R_c = inputVals.R_c;
else
    R_c = 900*10^3;
end
% Thickness of stagnant lid [m]
if isfield( inputVals , 'D_l' )
    D_l = inputVals.D_l;
else
    D_l = 200*10^3; % no idea what this value should properly be
end
% Silicate mantle density [kg / m^3]
if isfield( inputVals , 'rho_m' )
    rho_m = inputVals.R_m;
else
    rho_m = 3300; % cites Sohl et al. 2002
end
% Core density [kg / m^3]
if isfield( inputVals , 'rho_c' )
    rho_c = inputVals.rho_c;
else
    rho_c = 5000;
end
% Mantle thermal expansivity [1 / K]
if isfield( inputVals , 'alpha_m' )
    alpha_m = inputVals.alpha_m;
else
    alpha_m = 2*10^(-5); % cites Morschhauser et al. 2011
end
% Mantle heat capacity [J / K*kg]
if isfield( inputVals , 'c_pm' )
    c_pm = inputVals.c_pm;
else
    c_pm = 1150; % cites Hauck et al. 2006
end
% Mantle thermal diffusivity [m^2 / s]
if isfield( inputVals , 'kappa_m' )
    kappa_m = inputVals.kappa_m;
else
    kappa_m = 10^(-6); % cites Hauck et al. 2006
end
% Mantle viscocity [Pa * s]
if isfield( inputVals , 'nu' )
    nu = inputVals.nu;
else
    nu = 10^20.5; % geometric mean from Table 1
end
    
% gravitational constant [m^3 kg / s^2]
G = 6.67430*10^(-11);
% gravitational acceleration of the mantle
g_m = 4/3*pi*G*R_m*(rho_m + (rho_c-rho_m)*(R_c/R_m)^3); % cites Sohl et al. 2009a
% gravitational acceleration at CMB
g_c = 4/3*pi*G*rho_c*R_c;
% radius at bottom of stagnant lid
R_l = R_m - D_l;

%%% Calculating T_l
% constants
Theta = 2.9; % cites Morschhauser et al. 2011
R = 8.3144;

% temperature at the base of the stagnant lid
T_l = T_m - Theta*R*(T_m)^2 / E_m;

%%% Calculating T_b (1st pass, used for delta_i calculations)
% height between lower and upper boundary layers
%DeltaR = R_m - R_l - delta_t - delta_b;
DeltaR = R_l - R_c;% - delta_t - delta_b;
% deltas omitted due to circular dependence
% replaced R_m - R_l with R_l - R_c since it makes physical sense
% -Artyom

% temperature at the top  of the lower boundary layer
T_b = T_m + (alpha_m*g_m*T_m*DeltaR)/c_pm;

%disp(strcat('T_b = ',sprintf('%0.5e',T_b)))

for i = 1:5
    %%% Calculating delta_t
    % critical Rayleigh number
    Ra_crit = 450; % cites Choblet and Sotin 2000
    % delta T
    DeltaT = T_m - T_l + T_c - T_b;
    % Rayleigh number
    Ra = alpha_m * rho_m * g_m * DeltaT * (R_l - R_c)^3 / (kappa_m * nu); % B10\

    % thickness of upper boundary layer
    delta_t = (R_l - R_c) * (Ra_crit / Ra)^(1/3); % B9

    %%% Calculating delta_b
    % viscocity at average temperature in lower thermal boundary layer
    nu_b = nu*(T_c+T_b)/2; % cites Richter 1978
    % internal Rayleigh number
    Ra_i = alpha_m*rho_m*g_m*(T_m-T_int+T_c-T_b)*(R_m-R_c)^3/(kappa_m*nu); % cites Deschamps and Sotin 2000
    % Rayleigh number of lower thermal boundary layer
    Ra_bcrit = 0.28*Ra_i^(0.21); % cites Deschamps and Sotin 2000

    % thickness of lower boundary layer
    delta_b = ( (kappa_m * nu_b * Ra_bcrit) / (alpha_m * rho_m * g_c * (T_c - T_b)) )^(1/3); %B12

    %%% Calculating T_b
    % height between lower and upper boundary layers
    %DeltaR = R_m - R_l - delta_t - delta_b;
    DeltaR = R_l - R_c - delta_t - delta_b;
    % replaced R_m - R_l with R_l - R_c since it makes physical sense
    % -Artyom

    % temperature at the top  of the lower boundary layer
    T_b = T_m + (alpha_m*g_m*T_m*DeltaR)/c_pm;
    
    %disp(strcat('delta_t = ',sprintf('%0.5e',delta_t)))
    %disp(strcat('delta_b = ',sprintf('%0.5e',delta_b)))
    %disp(strcat('T_b = ',sprintf('%0.5e',T_b)))
end

%%% Calculating heat fluxes out of the mantle and the core
F_l = k_m * (T_m - T_l) / delta_t;
F_c = k_m * (T_c - T_b) / delta_b;