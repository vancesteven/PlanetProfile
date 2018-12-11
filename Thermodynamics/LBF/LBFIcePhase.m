function phase = LBFIcePhase(P,T,wt,salt)
% phase = getIcePhaseMgSO4(P,T,w)
% P in MPa
% T in K
% w in Wt%
% salt as |'NaCl'|, coming soon: |'MgCl2'|'MgSO4'|'Na2SO4'| 
%w in wt%


MW_H2O = 1e-3*18.01528; % in kg/mol
% iceII = load('iceII_IAPWS95_sp_G_fPT'); % coming soon
iceIII = load('iceIII_IAPWS95_sp_G_fPT'); 
iceV = load('iceV_IAPWS95_sp_G_fPT');
iceVI = load('iceVI_IAPWS95_sp_G_fPT');

out = GibbsIceIh({P T});                GIh = out.G; GIh(~GIh)=Inf;
                                        GII = Inf;
out = fnGval(iceIII.sp_G_fPT,[P T]);    GIII = out.G; GIII(~GIII)=Inf;
                                        GIV = Inf;
out = fnGval(iceV.sp_G_fPT,[P T]);      GV = out.G; GV(~GV)=Inf;
out = fnGval(iceVI.sp_G_fPT,[P T]);     GVI = out.G; GVI(~GVI)=Inf;

switch(salt)
    case 'MgCl2'
        error('MgCL2 is currently not supported')
    case 'MgSO4'
        error('MgSO4 is currently not supported')
    case 'NaCl'
        fluid = load('NaCl_LBF'); MW = 58.44;
        mm = WtPercent2Molality(wt,MW);        
        Gsalt = fnGval(fluid.sp_NaCl_8GPa,[P T mm],1e-3*MW);
        if isinf(Gsalt.rho)
            error('fluid P|T|m appear to be out of the bounds of the EOS')
        end
    case 'Na2SO4'
        error('Na2SO4 is currently not supported')
end

dmu = [MW_H2O*[GIh GII GIII GIV GV GVI] Gsalt.muw];
phase = find(dmu == min(dmu));

phase(phase==7)=0;
