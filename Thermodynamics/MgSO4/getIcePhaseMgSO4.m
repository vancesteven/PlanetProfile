function phase = getIcePhaseMgSO4(P,T,w)
% phase = getIcePhaseMgSO4(P,T,w)
% P in MPa
% T in K
% w in Wt%

%w in wt%
xH2O = 1-1./(1+(1-0.01.*w)./(0.01.*w).*120.3686./18.0142);

dmu = deltaMuILiqMgSO4(P,T,xH2O);

phase = find(dmu == min(dmu));

if dmu(phase) > 0
    phase = 0;
end

phase(phase==5) = 6;
phase(phase==4) = 5;

function dmu = deltaMuILiqMgSO4(P,T,x)
R = 8.314;
dmu = MuIceLiqH2O(P,T) - R./0.018.*T.*log(activity(P,T,x,-1.8e6,150,1.45e-4,-12,246)*x);% this line contains the coefficients from Vance et al. 2014 for MgSO4

function y = activity(P,T,x,w0,w1,w2,w3,To)
R = 8.314;
y = exp(((1-x).^2).*w0.*(1+w1.*tanh(w2.*P)).*(1+w3./(T-To).^2)./(R.*T));
