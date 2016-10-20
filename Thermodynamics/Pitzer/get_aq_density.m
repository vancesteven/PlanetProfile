function [rho_soln_gmL,rho_water_gmL,Vo,Vex] = get_aq_density(m_molal,P_bar,T_K,str_EOS)

W_salt.Na = 22.99;
W_salt.Mg = 24.31;
W_salt.Ca = 40.08;
W_salt.Cl = 35.4527;
W_salt.SO4 = (32.06+4*16.00);

VoType = 'Millero';
str_ion = fieldnames(m_molal);

mW = get_saltData(m_molal);
n_ions = length(str_ion);
lm = length(m_molal); 
lP = length(P_bar);
lT = length(T_K);

switch str_EOS
    case 'IAPWS'
        bar2GPa = 1e-4;
        rho_water_gmL = get_density(P_bar*bar2GPa,T_K);
    case 'Vance2008'
        bar2MPa = 0.1;
        for ij = 1:length(P_bar)
            for jk = 1:length(T_K)
                rho_water_gmL(ij,jk) = 1e-3*getWaterDensityVance2008(P_bar(ij)*bar2MPa,T_K(jk)-273.15);
            end
        end
    case 'Millero'
        for ij = 1:length(P_bar)
            V_water(ij,:) = getVwater_marion(P_bar(ij),T_K);
            rho_water_gmL(ij,:) = 1000*V_water(ij,:).^-1;
        end
end
Vw_m3 = 1000./rho_water_gmL;

mVo = zeros(lm,lP,lT);
Vo = mVo;
rho_soln_gmL = ones(lm,lP,lT);
Vex = rho_soln_gmL;
V_aq_m3 = rho_soln_gmL;

for hi = 1:lm
    mW = 0;
    for hii = 1:n_ions
         this_m = getfield(m_molal(hi),str_ion{hii});
         this_W = getfield(W_salt,str_ion{hii});
         mW = mW + this_m*this_W;
        switch VoType
            case 'supcrt92'
                mVo(hi,:,:) = squeeze(mVo(hi,:,:))+this_m*get_Xo_supcrt_fromFile(P_bar,T_K-273.15,str_ion{hii},'Vo');
            case 'Millero'
                Vo(hi,:,:) = squeeze(Vo(hi,:,:))+(getVo_Millero(P_bar',T_K,str_ion{hii}));
                mVo(hi,:,:) = squeeze(mVo(hi,:,:))+this_m*(getVo_Millero(P_bar',T_K,str_ion{hii}));
        end
    end
%     m_vec_molal = get_m_vec_molal(m_molal(hi));
    this_Vex = get_Vexs(m_molal(hi),P_bar,T_K); % this once used m_vec_molal
%     if length(T_K) == 1
%         this_V_aq_m3 = Vw_m3+mVo(hi,:,:)'+this_Vex';
%     else
%         this_V_aq_m3 = get_Vtotal(P_bar,T_K,Vw_m3,mVo(hi,:,:),this_Vex);
        this_V_aq_m3 = squeeze(Vw_m3)+squeeze(mVo(hi,:,:))+squeeze(this_Vex);
%     end
    % solution density calculation
    rho_soln_gmL(hi,:,:) = (1000+mW)./this_V_aq_m3; % rho in gm/mL
    Vex(hi,:,:) = squeeze(this_Vex);
    V_aq_m3(hi,:,:) = this_V_aq_m3;
end