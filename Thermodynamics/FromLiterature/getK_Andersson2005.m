function K_W_mK = getK_Andersson2005(P_MPa,T_K,varstr,typestr)

if typestr=='T'
    K_W_mK = kT(T_K,varstr);
elseif typestr=='P'
    K_W_mK = kP(P_MPa,varstr);
end


function K_W_mK = kT(T_K,varstr)
% temperature dependence (table 1)
% Phase D/W m1/K11x X p/MPa Tmin Tmax /K Inaccuracy (%)

Ih = [630 0.995 0.1 40 180 5];
Ic =[ 192 0.813 0.1 150 200 6];
II =[ 695 1.097 240 120 240 4];
III =[ 93.2 0.822 240 180 250 4];
% IV ? ? ? ? ?
V =[ 38.0 0.612 530 240 270 4];
VI=[ 50.9 0.612 1000 135 250 4];
VII=[ 320 0.821 2400 275 300 4];
VIII=[ 15 700 1.417 2400 240 270 4];
IX=[ 164 0.879 240 120 160 4];
% X ? ? ? ? ?
XI=[ 994 1.041 0.1 50 72 6];
% Metastable ice
% (ice XII)17
XII=[80.1 0.77 100 115 140 4];
% ASW ? ? ? ? ?
LDA=[ 17.5 0.55 0.1 80 130 4];
HDA=[ 0.528 0.026 0.1 75 120 4];
VHDA=[ 0.257 0.20 1000 145 155 6];
if isvarname(varstr)
    eval(['coeffs = ' varstr ';']);
else
    coeffs = [0 0];
end
Dwm = coeffs(1);
xT = coeffs(2);
K_W_mK = Dwm*T_K.^-xT;


function K_W_mK = kP(P_MPa,varstr)
% pressure dependence (table 3)
% Phase E F/GPa1 pmin pmax/GPa T/K
Ih=[ 1.60 -0.44   0 0.50 130];
Ic=[ 1.28 -0.34   0 0.50 130];
II=[ 1.25 0.2 0   0.24 120];
III=[-0.02 0.2  0.2 0.35 240];
% IV ? ? ? ?
V=[ 0.16 0.2      0.35 0.60 246];
VI=[ 0.37 0.16    0.7 2.0 246];
VII=[ 0.65 0.2    2.0 2.4 286];
VIII=[ 1.38 0.2   2.0 2.4 246];
% IX ? ? ? ?
% X ? ? ? ?
XI=[ 2.67 0.24    0.08 0.16 58];
% Metastable ice (ice XII)17 
XII=[0.73 0.09      0.1 0.8 115];
% ASW ? ? ? ?
LDA=[0.17 -0.63   0 0.35 130];
HDA130K18=[-0.52 0.19 0 0.50 130];
% VHDA ? ? ? ?
if isvarname(varstr)
    eval(['coeffs = ' varstr ';']);
end
Ek = coeffs(1);
F_GPam1 = coeffs(2);

%p dependence
K_W_mK = Ek + F_Gpam1*P_MPa*1e-3;

