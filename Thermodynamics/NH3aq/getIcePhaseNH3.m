function phase = getIcePhaseNH3(P,T,w)
% phase = getIcePhaseNH3(P,T,w)
% P in MPa
% T in K
% w in Wt%

% mNH3 = 14.0067 + 3*1.00794 = 17.03052

%w in wt%
xH2O = 1-1./(1+(1-0.01.*w)./(0.01.*w).*17.03052./18.0142); % mole fraction of water

dmu = deltaMuILiqNH3(P,T,xH2O);

phase = find(dmu == min(dmu));

if dmu(phase) > 0
    phase = 0;
end

phase(phase==5) = 6;
phase(phase==4) = 5;

function dmu = deltaMuILiqNH3(P,T,x)
R = 8.314;
dmu = MuIceLiqH2O(P,T) - R./0.018.*T.*log(activity(P,T,x,-676.9e3,0.4,6.0e-3,4.62e-4,-113.0)*x);% this line contains the coefficients from Choukroun and Grasset 2010

function y = activity(P,T,x,w0,w1,w2,w3,w4)
R = 8.314./0.018;
y = exp(((1-x).^2).*w0.*(1+w1.*tanh(w2.*P)).*(1+w3.*T+w4./T)./(R.*T));