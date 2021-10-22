function n = WtPercent2Molality(wtpcnt,W_salt)
%  n = WtPercent2Molality(wtpcnt,W_salt)
% n = 1000*wtpcnt/(100-wtpcnt)/W_salt;
% WMgSO4 = 120.4 g mol-1
%WNH3 = 17.03026; % g mol-1 
Wwater = 18.015268; % g mol-1 

n = 1000*wtpcnt./(100-wtpcnt)./W_salt;