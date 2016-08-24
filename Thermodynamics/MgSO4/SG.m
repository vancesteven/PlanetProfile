function [Tm] = SG(i,Pm)
%%% Simon-Glatzel equations for ices Ih-III-V-VI from Choukroun and Grasset (2007)
%%% Tm (K), melting temperature
%%% Pm (MPa), melting pressure

if not(ismember(i,[1 3 5 6]))
    error('SG equation error - parameters available for ices Ih, III, V and VI only')
end
    
%%% Simon-Glatzel equations parameters
%%% i(phase): 1(Ih), 3(III), 5(V), 6(VI)
P0 = [6.11657e-4   0   209.50   0   355.00   618.40] ;
T0 = [273.15       0   251.15   0   256.43   272.73] ;
a  = [-414.5       0   101.10   0   373.60   661.40] ;
c  = [8.3800       0   42.860   0   8.6600   4.6900] ;

% Pm = P0(i) + a(i) * (((Tm/T0(i))^c(i))-1) ;

Tm = T0(i) * ((((Pm-P0(i))/a(i))+1)^(1/c(i))) ;

end