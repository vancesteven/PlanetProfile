function [Vo,Cp_J_kg_K] = getVo_Millero(P,T,str_ion)
% T in K
% P in bars
%adopted from Marion (Geochim Cosmochim Acta, 2005), Table 1 - fits to Millero (Chem Ocean (8), 1983)
%Vo = Vo_T - P*beta(T);

[Vo_T,beta] = getVo_TbetaMillero(T,str_ion);
Vo = ones(length(P),1)*Vo_T - (P-1)*beta; % marion 2005 equation 24
%==========================================================================
function [Vo_T,beta] = getVo_TbetaMillero(T,str_ion)
switch str_ion
    case 'Na'
        aV = [-7.68e-4 5.2876e-1 -9.0589e1];
        ab = [-1.56e-7 1.36923e-4 -3.0891e-2];
    case 'K'
        aV = [7.36e-4 4.9128e-1 7.2019e1];
        ab = [5.6e-8 7.613e-6 -1.0616e-2];
    case 'Fe'
        aV =[1.600e-3 9.6088e-1 1.65407e2];
        ab = [8.44e-7 -4.8948e-4 6.2899e-2];
    case 'Mg'
        aV = [-1.600e-3 9.6088e-1 -1.65407e2];
        ab = [8.44e-7 -4.8948e-4 6.2899e-2];
    case 'Ca'
        aV = [-1.256e-3 7.9195e-1 -1.42301e2];
        ab = [8.44e-7 -4.8948e-4 6.2899e-2];
    case 'Cl'
        aV = [-1.264e-3 7.8012e-1 -1.02412e2];
        ab = [-1.644e-6 1.0101e-3 -1.5576e-1];
    case 'OH'
        aV = [-3.536e-3 2.16412e0  -3.34884e2];
        ab = [-1.160e-6 7.6471e-4 -1.2874e-1];
    case 'HSO4'
        aV = [-2.248e-3 1.41308e0 -1.85898e2];    
        ab = [-1.160e-6 7.6471e-4 -1.2737e-1];
    case 'NH4'
        aV = [1.058741e-3   -5.920729e-1    1.002813e-2];
        ab = [-4.895105e-8  -1.181049e-5    9.069507e-3]; % Marion et al. 2012, valid from 273-323
    case 'NO3'
        aV = [-1.984e-3 1.27266e0 -1.74089e2];
        ab = [-6.08e-7 4.2295e-4  -7.1475e-2];
    case 'HCO3'
        aV = [-2.248e-3 1.41308e0 -1.97188e2];
        ab = [-1.160e-6 7.6471e-4 -1.2763e-1];
    case 'CO3'
        aV = [-4.064e-3 2.52016e0 -3.93904e2];
    	ab = [-3.056e-6 1.8935e-3 -3.0024e-1];
    case 'SO4'
        aV = [-3.808e-3 2.3643    -3.52433e2];
        ab = [-3.056e-6 1.8935e-3 -2.9963e-1];
    case 'CO2' 
        aV = [0 0 33.6];
    case 'O2' 
        aV = [0 0 32.0];
    case 'CH4'
        aV = [0 0 37.4];        
end
Vo_T = polyval(aV,T);
beta = polyval(ab,T)*1/2; % /2 gives a better fit to density in pressure, but only in the T  = 0:50 oC range
 
%==========================================================================
function Cp_J_kg_K = get_CpMillero2005(P,T_K)
coeffs = [-0.00362131
           0.0001355158
           -2.79789e-6
            3.611949e-08
           -2.840693e-10
            1.342793e-12
            -3.48798e-15
            3.82206e-18];
cp0_J_g_K= 4.21713 +polyval(coeffs(end:-1:1),T_K);

        