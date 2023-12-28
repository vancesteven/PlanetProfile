function phase = getIcePhaseMgSO4(P,T,w)
% phase = getIcePhaseMgSO4(P,T,w)
% P in MPa
% T in K
% w in Wt%

%w in wt%
xH2O = 1-1./(1+(1-0.01.*w)./(0.01.*w).*120.3686./18.0142);

dmu = deltaMuILiqMgSO4(P,T,xH2O);

nTs = length(T);
phase = zeros(1,nTs);
for iT=1:nTs
    phase(iT) = find(dmu(:,iT) == min(dmu(:,iT)));
    
    if all(dmu(:,iT) > 0)
        phase(iT) = 0;
    end
end

phase(phase==5) = 6;
phase(phase==4) = 5;
phase = unit8(phase); % Cast to small integer

function dmu = deltaMuILiqMgSO4(P,T,x)
Rsp = 8.314/0.018;
[dmuIce, ~, ~] = MuIceLiqH2O(P,T);
dmu = zeros(5,length(T));
for iT=1:length(T)
    dmu(:,iT) = dmuIce(:,:,iT) - (activity(P,T(iT),x,-1.8e6,150,1.45e-4,-12,246) + Rsp.*T.*log(x));% this line contains the coefficients from Vance et al. 2014 for MgSO4
end

function y = activity(P,T,x,w0,w1,w2,w3,To)
y = ((1-x).^2).*w0.*(1+w1.*tanh(w2.*P)).*(1+w3./(T-To).^2);