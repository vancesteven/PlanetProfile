function Cp_J_kg_K = CpH2O_Choukroun(P_MPa,T_K,ind)
% from Choukroun and Grasset 2010
% inds = getInds(P_MPa,T_K,ind);
% addition for methane clathrates from Ning et al 2014;
Cp_J_kg_K = getCp2010(T_K,ind);

function Cp_J_kg_K = getCp2010(T_K,inds)
%  Ice Ih, Ice II, Ice III, Ice V, Ice VI, methane clath
c0 = [4190 74.11 2200 820  700 940 2150]; 
c1 = [9 7.56 0 7  7.56 5.5 3.19];

Cp_J_kg_K = polyval([c1(inds) c0(inds)],T_K);