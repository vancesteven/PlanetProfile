function T_l_K = SulfurRichLiquidusRuckrieman(p_GPa,x_s)
% Calculation of liquidus based on Ruckrieman et al. 2018 (A2)
% p_GPa = pressure
% x_s = sulfut wt %

p = p_GPa*1e9;

T_eut = 1265-11.15e-9*(p-3e-9);
x_eut = 0.11 + 0.187*exp(-0.065e-9*p);
T_FeS = 1500.5+29.6e-9*p - 0.247e-18*p^2;
T_l_K = T_FeS + (T_FeS-T_eut)./(0.3647-x_eut).*(x_s-0.3647);