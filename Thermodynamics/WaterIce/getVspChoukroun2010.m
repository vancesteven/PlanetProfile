function V_cm3_kg = getVspChoukroun2010(P_MPa,T_K,inds)

% liquid, Ice Ih, Ice II, Ice III, Ice V, Ice VI
Tref_K =    [400    273.16   238.45     256.43      273.31      356.15];
Vo_dm3_kg = [0.815   1.086   0.8425     0.855       0.783       0.743];

a0 = [0.100     0.019       0.060       0.0375      0.005       0.024];
a1 = [0.0050    0.0075      0.0070      0.0203      0.0100      0.002];
b0 = [1.000     0.974       0.976       0.951       0.977       0.969];
b1 = [0.2840    0.0302      0.0425      0.097       0.1200      0.05];
%modifications following CG_unknownSH from Bollengier -- no correct
%pressure is found for the ice-water interface with this line:
%b2 = [0.00136    0.00395     0.0022      0.00102     0.0016      0.0002];
b2 = [0.00136   0.00395     0.0022      0.00200     0.0016      0.00102];

Tdiff_K = T_K-Tref_K(inds);
EpsT = 1 + a0(inds).*tanh(a1(inds).*Tdiff_K);
EpsP = b0(inds) + b1(inds).*(1-tanh(b2(inds).*P_MPa));
V_cm3_kg = Vo_dm3_kg(inds).*EpsT.*EpsP;
