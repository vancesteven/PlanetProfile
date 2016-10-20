function Vex2 = get_Vex2Pitz(m_molal,P_bar,T_K)
% Following Marion et al (2005).  T in K, P in bars
% Calculate Vex =
% Av*I/b*ln(1+b*I^0.5)+2*R*T*sum(sum(mc*ma*(Bcav+sum(mc*zc)*Ccav)))
% Krumgalz (1995) has two additional terms in the last summation.  Marion
% (2005) truncates them

% Bcav = Bcav0+Bcav1*g(alpha1*I^0.5)_Bcav2*g(alpha2*I^0.5)
% Bcav' = Bcav1*dg(x1)/dI+Bcav2*dg(x2)/dI
% Polynomial coefficients from Marion et al (2005) Table 3
lm = length(m_molal);
Vex2 = zeros(lm,length(P_bar),length(T_K));

for hi = 1:lm
    if isfield(m_molal,'SO4')
        if isfield(m_molal,'Na')
            z_cat = 1;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'Na','SO4');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).Na,m_molal(hi).SO4,P_bar,T_K,B,C,z_cat);
        end
        if isfield(m_molal,'Mg')
            z_cat = 2;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'Mg','SO4');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).Mg,m_molal(hi).SO4,P_bar,T_K,B,C,z_cat);
        end
        if isfield(m_molal,'Fe')
            z_cat = 2;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'Fe','SO4');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).Fe,m_molal(hi).SO4,P_bar,T_K,B,C,z_cat);
        end
    end
    if isfield(m_molal,'Cl')
        if isfield(m_molal,'Na')
            z_cat = 1;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'Na','Cl');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).Na,m_molal(hi).Cl,P_bar,T_K,B,C,z_cat);
        end
        if isfield(m_molal,'Mg')
            z_cat = 2;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'Mg','Cl');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).Mg,m_molal(hi).Cl,P_bar,T_K,B,C,z_cat);
        end
        if isfield(m_molal,'K')
            z_cat = 1;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'K','Cl');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).K,m_molal(hi).Cl,P_bar,T_K,B,C,z_cat);
        end
        if isfield(m_molal,'Fe')
            z_cat = 2;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'Fe','Cl');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).Fe,m_molal(hi).Cl,P_bar,T_K,B,C,z_cat);
        end
        if isfield(m_molal,'Ca')
            z_cat = 2;
            [B C] = get_PitzerCoeffs(m_molal(hi),P_bar,T_K,'Ca','Cl');
            Vex2(hi,:,:) = squeeze(Vex2(hi,:,:)) + get_Vex2(m_molal(hi).Ca,m_molal(hi).Cl,P_bar,T_K,B,C,z_cat);
        end
    end
end
VexNaNs =  find(isnan(Vex2));
Vex2(VexNaNs) = 0;
%=================================================m=========================
function Vex2 = get_Vex2(m_cat,m_an,P_bar,T_K,B,C,z_cat)
R = 83.1451; % cm3*bar/mol/K
lP = length(P_bar);  lT = length(T_K);
for ij = 1:lP
    for jk = 1:lT 
       Vex2(ij,jk) = 2*R*T_K(jk).*m_cat.*m_an.*(B(:,jk)'+m_cat*z_cat*C(jk));
    end
end
%=================================================m=========================
function [B C] = get_PitzerCoeffs(m_molal,P_bar,T_K,str_cat,str_an)
% from Krumgalz (2000), page 1125, alphax_n_m || (kg/mol)^0.5
alpha1_1x = 2.0;
alpha1_22 = 1.4;
alpha2_22 = 12; % this one is only valid for 2:2 electrolytes (like MgSO4), since Bcav2 is assumed zero for 1:x and x:1 electrolytes

lT = length(T_K);
switch str_an
    case 'SO4'
        switch str_cat
            case 'Na'            
                z_cat_sq = 1; z_an_sq = 4;
                I = 0.5*(m_molal.Na*z_cat_sq+m_molal.SO4*z_an_sq);
                B0 = polyval([-1.698e-6 5.589854e-4],T_K);
                B1 = polyval([3.717169e-7 -2.218483e-4 3.323357e-2],T_K);
                B2 = zeros(lT);
                C = polyval([2.106e-7 -6.561572e-5],T_K);
                gaI2 = 0;
				gaI1 = 2/alpha1_1x^2./I.*(1-(1+alpha1_1x*I.^0.5).*exp(-alpha1_1x.*I.^0.5));
            case 'Mg'
                z_cat_sq = 4; z_an_sq = 4;
                I = 0.5*(m_molal.Mg*z_cat_sq+m_molal.SO4*z_an_sq);
                B0 = polyval([-4.330e-7 1.789039e-4],T_K);
                B1 = polyval([3.358883e-7 -2.055400e-4 3.157137e-2],T_K);
                B2 = polyval([-2.570e-4 9.096288e-2],T_K);
                C = polyval([-2.629e-8 8.232980e-6],T_K);
                gaI2 = 2/alpha2_22^2./I.*(1-(1+alpha2_22*I.^0.5).*exp(-alpha2_22*I.^.5));
				gaI1 = 2/alpha1_22^2./I.*(1-(1+alpha1_22*I.^0.5).*exp(-alpha1_22.*I.^0.5));
            case 'NH4'
                z_cat_sq = 1; z_an_sq = 4;
                I = 0.5*(m_molal.NH4*z_cat_sq+m_molal.SO4*z_an_sq);
                B0 = polyval([8.378e-4  -2.107e-1],T_K);
                B1 = polyval([2.264e-2  -6.089e0],T_K);
                B2 = zeros(lT);
                C = polyval([-9.424e-5  2.768e-2],T_K);
                gaI2 = 0;
				gaI1 = 2/alpha1_1x^2./I.*(1-(1+alpha1_1x*I.^0.5).*exp(-alpha1_1x.*I.^0.5));
            case 'Fe'
                z_cat_sq = 4; z_an_sq = 4;
                I = 0.5*(m_molal.Fe*z_cat_sq+m_molal.SO4*z_an_sq);
                B0 = polyval([1.196667e-6 5.246895e-4],T_K);
                B1 = polyval([5.846667e-6 3.368357e-3],T_K);
                B2 = zeros(lT);
                C = zeros(lT);
                gaI2 = 2/alpha2_22^2./I.*(1-(1+alpha2_22*I.^0.5).*exp(-alpha2_22*I.^.5));
				gaI1 = 2/alpha1_22^2./I.*(1-(1+alpha1_22*I.^0.5).*exp(-alpha1_22.*I.^0.5));
        end
    case 'Cl'
        switch str_cat
            case'Na'
                z_cat_sq = 1; z_an_sq = 1;
                I = 0.5*(m_molal.Na*z_cat_sq+m_molal.Cl*z_an_sq);
                B0 = polyval([-3.2412e-7 1.088468e-4],T_K);
                B1 = polyval([6.911453e-8 -4.135834e-5 6.193806e-3],T_K);
                B2 = zeros(lT);
                C = polyval([1.514940e-8 -5.174534e-6],T_K);
            case 'Mg'
                z_cat_sq = 4; z_an_sq = 1;
                I = 0.5*(m_molal.Mg*z_cat_sq+m_molal.Cl*z_an_sq);
                B0 = polyval([6.920e-9 -4.168596e-6 6.446574e-4],T_K);
                B1 = polyval([-5.618320e-6 1.624497e-3],T_K);
                B2 = zeros(lT);
                C = -5.567e-7*T_K;
            case 'K'
                z_cat_sq = 1; z_an_sq = 1;
                I = 0.5*(m_molal.K*z_cat_sq+m_molal.Cl*z_an_sq);
                B0 = polyval([1.782e-7 6.593033e-5],T_K); 
                B1 = polyval([3.224932e-5 5.215803e-8 4.987530e-3],T_K);
                B2 = zeros(lT);
                C = 7.112e-7*T_K;
            case 'Fe'
                z_cat_sq = 4; z_an_sq = 1;
                I = 0.5*(m_molal.Fe*z_cat_sq+m_molal.Cl*z_an_sq);
                B0 = polyval([2.42e-7 1.19712e-4],T_K);
                B1 = polyval([1.80e-6 2.40397e-3],T_K);
                B2 = zeros(lT);
                C = zeros(lT);
            case 'Ca'
                z_cat_sq = 4; z_an_sq = 1;
                I = 0.5*(m_molal.Ca*z_cat_sq+m_molal.Cl*z_an_sq);
                B0 = polyval([-4.293e-7 1.409126e-4],T_K);
                B1 = polyval([2.612103e-6 -7.785779e-4 7.702445e-2],T_K);
                B2 = zeros(lT);
                C =  polyval([2.216e-8 -6.723184e-6],T_K);
        end
            % All Cl are 1x
                gaI2 = 0;
				gaI1 = 2/alpha1_1x^2./I.*(1-(1+alpha1_1x*I.^0.5).*exp(-alpha1_1x.*I.^0.5));
end
for jk = 1:lT
    if ~isempty(gaI1)
        B(:,jk) = B0(jk)+B1(jk).*gaI1+B2(jk).*gaI2;
    else
        B(:,jk) = 0;
    end
end