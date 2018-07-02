function [V]=iceVI_PT_EOS_bezacier(P,T)
%Ice VI BM2 PVT-EOS Bezacier et al. (2014) (best > 300K) (Implemented by B. Journaux)
% Rerturns pressure in GPa for T (K) and V (cm^3/mol)

% Ice VI BM2 PVT-EOS Bezacier et al. (2014)
K0 = 14.05168; % GPa
Kp = 4;
V0 = 14.17062479 ; % cm^3 / mol
alpha = 14.63878 * 10.^(-5); % K^-1
T0 = 300;

%% BM2 calculation :
% Thermal volume :
V0_T = V0.* exp(alpha.*(T-T0)); % V0(T) calculation from V(P0,T0) and alpha
% Converge
fun = @(x)P - 3* K0.*((V0_T./x).^(2./3)-1)/2.* (1 + 2.* ((V0_T/x).^(2./3)-1)/2.).^(5/2).* (1 + 3/2.* (Kp-4).* ((V0_T./x).^(2./3)-1)./2.); 
x0 = 14;
options = optimoptions('fsolve','Display','none');
V = fsolve(fun,x0,options);