function wtpcnt = Molality2WtPercent(molality,W_salt)
% Converts a given molality for an aqueous species of molar mass W_salt to weight percent 
% USAGE: wtpcnt = Molality2WtPercent(molality,W_salt)
% e.g., W_salt:
% W_Mg = 24.3050, W_S = 32.065, W_O = 15.9994 => W_MgSO4 = 120.3
% W_Na  = 22.98977, W_Cl = 35.453
% W_Na2SO4 = 142.04;
wtpcnt = 100*molality*W_salt./(1000+molality*W_salt);
