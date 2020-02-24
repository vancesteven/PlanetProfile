function PPGanymedeCMExploration
%PPGanymede
Planet.name='Ganymede';
Planet.rho_kgm3 = 1936; % Schubert et al. 2004: 1942.0±4.8 claimed to be 4x more accurate than Anderson 1996 of 1936 ± 22
Planet.R_m = 2634.1e3;
Planet.M_kg =1.4819e23;
Planet.gsurf_ms2 = 1.428; 
Planet.Tsurf_K = 110; 
Planet.Psurf_MPa = 0; 



%% Porosity of the rock
Planet.POROUS_ROCK = 0;


%% Seismic
Seismic.LOW_ICE_Q = 1; % divide Ice Q value by this number
Seismic.QScore = 1e4;
Seismic.SMOOTH_VROCK = 1; % smooth over N neighboring rows and columns in vp and vs

%Attenuation Parameters Based on those Described in Cammarano et al. 2006
% ice I
Seismic.B_aniso_iceI = 0.56;
Seismic.gamma_aniso_iceI = 0.2;
Seismic.g_aniso_iceI = 22; %C2006
% ice II
Seismic.B_aniso_iceIII = 0.56;
Seismic.gamma_aniso_iceIII = 0.2;
Seismic.g_aniso_iceIII = 30; 
% ice III
Seismic.B_aniso_iceIII = 0.56;
Seismic.gamma_aniso_iceIII = 0.2;
Seismic.g_aniso_iceIII = 25; 
% ice V
Seismic.B_aniso_iceV = 0.56;
Seismic.gamma_aniso_iceI = 0.2;
Seismic.g_aniso_iceV = 27; 
% ice VI
Seismic.B_aniso_iceVI = 0.56;
Seismic.gamma_aniso_iceVI = 0.2;
Seismic.g_aniso_iceVI = 28; 
% mantle
Seismic.B_aniso_mantle = 0.56;
Seismic.gamma_aniso_mantle = 0.2;
Seismic.g_aniso_mantle = 30; %C2006

%% Model Parameters
Params.foursubplots =1;
Params.HOLD = 0; % overlay previous run
Params.Legend=0;
Params.Pseafloor_MPa = 2000;
Params.LegendPosition = 'southeast';
Params.ylim = [925 1375];
Params.nsteps_iceI = 50;
Params.nsteps_ocean = 600; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 100;
Params.nsteps_core = 10;
Params.savefigformat = 'epsc';
Params.wref=[0 5 10 15];
Params.BOTTOM_ICEIII=0;
Params.BOTTOM_ICEV=0;


Params.colororder = 'cbmkgrm';
Params.Temps = [250 252.5 255 260 265 270 273];

Params.wrefLine = '--';


Planet.Cuncertainty = 0.0028;
%note: Schubert et al. 2004 cite the Anderson C/MR2 as 0.3115±0.0028.  This
%is incorrect, as the value cited everywhere else is consistent with the
%Anderson et al. (1996) value of C/MR2=0.3105\pm0.0028 used here
Planet.rhoFe = 8000; %8000
Planet.rhoFeS = 5150; %5150

%  Interior constraints imposed in Vance et al. 2014
mSi = 28.0855; mS = 32.065; mFe = 55.845; mMg = 24.305;
xOl = 0.44; % percentage of olivine - Javoy (1995) - Fe/Si = 0.909 Mg/Si = 0.531, Mg# = 0.369
%mOl = 2*((1-0.369)*58.85+0.369*24.31)+28.0855+4*16=184.295
%mPx = 2*((1-0.369)*58.85+0.369*24.31+28.0855+3*16) =244.3805
xSi = (xOl+2*(1-xOl))*mSi/(xOl*184.295+(1-xOl)*244.3805); % mass fraction of sulfur in silicates
M_Earth_kg = 5.97e24;
xSiEarth = 0.1923; % Javoy in kg/kg in Exoplanets paper20052006-xSiSolar only in mole
xK = 1; %enrichment in K
Hrad0 = 24e12*xSi/xSiEarth/M_Earth_kg;


%% Mantle Heat
%cold case  
Planet.EQUIL_Q = 0;
Planet.kr_mantle = 4; % rock conductivity (Cammarano et al. 2006, Table 4)
Planet.Qmantle_Wm2 = 2.5e11/4/pi/Planet.R_m^2; 
Planet.QHmantle = 0;
%hot case Qm = 2.1e11+8.5e11; %W
% Qmantle = 1.3e11; 
% QHmantle = 8.5e11;



Params.LineStyle='--';



% these need to be adjusted.
Planet.xFeS_meteoritic = 0.0676; %CM2 mean from Jarosewich 1990
Planet.xFeS = 0.2; %0.25
Planet.xFe_core = 0.0463 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.rho_sil_withcore_kgm3 = 3730;

silicates = {
% % 'CI_hhph_DEW17_fluid_nomelt_685.tab'  
% 'CI_hhph_DEW17_nofluid_nomelt_685.tab'
% % 'CM_hhph_DEW17_fluid_nomelt_685.tab'  
% 'CM_hhph_DEW17_nofluid_nomelt_685.tab'
% % 'CV_hhph_DEW17_fluid_nomelt_685.tab'  
% 'CV_hhph_DEW17_nofluid_nomelt_685.tab'
% % 'H_hhph_DEW17_685_fluid_nomelt.tab'   
% 'H_hhph_DEW17_685_nofluid_nomelt.tab'   
% % 'L_hhph_DEW17_685_v1_fluid_nomelt.tab'  
% 'L_hhph_DEW17_685_v1_nofluid_nomelt.tab'
% 'Simple_CI_HS_green_PP.tab'
% 'Simple_CM_HS_green_PP.tab'
% 'Simple_CV_HS_green_PP.tab'
'CM_hydrous_differentiated_Ganymede_Core100Fe0S_excluding_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core95Fe5S_excluding_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core90Fe10S_excluding_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core85Fe15S_excluding_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core80Fe20S_excluding_fluid_properties.tab'
};
silicates_wFluids = {
'CM_hydrous_differentiated_Ganymede_Core100Fe0S_including_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core95Fe5S_including_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core90Fe10S_including_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core85Fe15S_including_fluid_properties.tab'
'CM_hydrous_differentiated_Ganymede_Core80Fe20S_including_fluid_properties.tab'    
};

mphases = {
    'CM_hydrous_differentiated_Ganymede_Core100Fe0S_wt_percent_phases.tab'
    'CM_hydrous_differentiated_Ganymede_Core95Fe5S_wt_percent_phases.tab'
    'CM_hydrous_differentiated_Ganymede_Core90Fe10S_wt_percent_phases.tab'
    'CM_hydrous_differentiated_Ganymede_Core85Fe15S_wt_percent_phases.tab'
    'CM_hydrous_differentiated_Ganymede_Core80Fe20S_wt_percent_phases.tab'    
    };
mfluids = {
    'CM_hydrous_differentiated_Ganymede_Core100Fe0S_mass_fraction_Fluid.tab'
    'CM_hydrous_differentiated_Ganymede_Core95Fe5S_mass_fraction_Fluid.tab'
    'CM_hydrous_differentiated_Ganymede_Core90Fe10S_mass_fraction_Fluid.tab'
    'CM_hydrous_differentiated_Ganymede_Core85Fe15S_mass_fraction_Fluid.tab'
    'CM_hydrous_differentiated_Ganymede_Core80Fe20S_mass_fraction_Fluid.tab'
};
mvolumes = {
    'CM_hydrous_differentiated_Ganymede_Core100Fe0S_vol_percent_phases.tab'
    'CM_hydrous_differentiated_Ganymede_Core80Fe20S_vol_percent_phases.tab' 
    'CM_hydrous_differentiated_Ganymede_Core85Fe15S_vol_percent_phases.tab'
    'CM_hydrous_differentiated_Ganymede_Core90Fe10S_vol_percent_phases.tab'
    'CM_hydrous_differentiated_Ganymede_Core95Fe5S_vol_percent_phases.tab'
};
cores={
  '100Fe0S_core.tab'  
  '95Fe5S_core.tab'  
  '90Fe10S_core.tab'  
  '85Fe15S_core.tab'  
  '80Fe20S_core.tab'    
 };

% [125 100 75 50]
Params.NOPLOTS = 1; %allows user to limit recreating plots & figures after each run

% Planet.Cmeasured = 0.3115;
cmr2 = 0.3115+[-0.0028 0 0.0028];
% cmr2 = 0.3115+[0.0028];
xS = [0 .05 0.10 0.15 0.20];
Planet.Cmeasured = 0.3115;
Planet.xFeS = 0.2; %0.25

Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
% rmfield(Planet.Ocean,'fnTfreeze_K');
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

Params.CALC_NEW =1; % Set CALC_NEW options to 0 to re-use profile data when possible. It is recommended to keep CALC_NEW=1 except when intermediate parameters such as layer thicknesses will not change between runs.
<<<<<<< HEAD
Params.CALC_NEW_REFPROFILES=0;
=======
Params.CALC_NEW_REFPROFILES=1;
>>>>>>> ee144b8c4a5cdb93c30a88553ad5b735300303c1
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY=1;
Planet.FeCore=true;

Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [254.75 259.05 263.1 266]; % 10 Wt% temperatures at the bottom of the Ice Ih
rhos = [3524 3522 3523 3523 3522
        3529 3528 3528 3526 3526
        3528 3532 3531 3529 3529];
for ic = 1:length(cmr2)
    Planet.Cmeasured = cmr2(ic);
%     for iS = 1:length(xS)
    for iS = length(xS)
        Planet.xS = xS(iS);
        Seismic.mantleEOS = silicates{iS};
        Seismic.coreEOS = cores{iS};
        Seismic.mantlefluidsEOS = silicates_wFluids{iS};
        Seismic.mfluids = mfluids{iS};
        Seismic.mphases = mphases{iS};
        Seismic.mvolumes = mvolumes{iS};
        disp(['input S%: ' num2str(Planet.xS)])
        disp(['silicates: ' silicates{iS} ])
        disp(['cores:     ' cores{iS} ])
        if rhos(ic,iS)
            disp([ic iS])
            Planet.rho_sil_withcore_kgm3 = rhos(ic,iS);
            PlanetProfile(Planet,Seismic,Params);
            Params.CALC_NEW =0; % only need to do this once. CALC_NEW only pertains to the volatile part. The rest of the calculation is quick and so is done every time.
        end
    end
end

<<<<<<< HEAD
Params.CALC_NEW =0;
Params.CALC_NEW_REFPROFILES=0;
=======


Params.CALC_NEW =1;
Params.CALC_NEW_REFPROFILES=1;
>>>>>>> ee144b8c4a5cdb93c30a88553ad5b735300303c1
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY=0;
Params.HOLD = 0; % overlay previous run
Planet.FeCore=true;
Planet.Ocean.w_ocean_pct=0;  Planet.Tb_K = [256.5  260.6  264.4 267.7]; % pure water, temperatures at the bottom of the Ice Ih
Params.LineStyle='-';
% PlanetProfile(Planet,Seismic,Params)
rhos = [3529 3529 3528 3527 3525
        3524 3531 3529 3528 3527
        3525 3533 3531 3529 3528];
for ic = 1:length(cmr2)
    Planet.Cmeasured = cmr2(ic);
    for iS = 1:length(xS)
        Planet.xS = xS(iS);
        Seismic.mantleEOS = silicates{iS};
        Seismic.coreEOS = cores{iS};
        if rhos(ic,iS)
            disp([ic iS])
            Planet.rho_sil_withcore_kgm3 = rhos(ic,iS);
            PlanetProfile(Planet,Seismic,Params);
            Params.CALC_NEW =0; % only need to do this once. CALC_NEW only pertains to the volatile part. The rest of the calculation is quick and so is done every time.
        end
    end
end

<<<<<<< HEAD
Params.CALC_NEW =0;
Params.CALC_NEW_REFPROFILES=0;
=======

Params.CALC_NEW =1;
Params.CALC_NEW_REFPROFILES=1;
>>>>>>> ee144b8c4a5cdb93c30a88553ad5b735300303c1
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY=0;
Params.HOLD = 1; % overlay previous run
Planet.FeCore=true;
Params.LineStyle='-';
% Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab'; % this needs to be substituted with the version that conserves Planet.rho_sil_withcore_kgm3 = 3518;
Seismic.mantleEOS = 'CI_anhydrous_differentiated_Ganymede_including_fluid_properties.tab';Planet.rho_sil_withcore_kgm3 = 3370;
Planet.Ocean.w_ocean_pct=1;  Planet.Tb_K = [256.5  260.6  264.4 267.7]; % pure water, temperatures at the bottom of the Ice Ih
% PlanetProfile(Planet,Seismic,Params)
rhos = [3352 3355 3353 3360 3363
        3330 3333 3330 3338 3340
        3308 3309 3307 3313 3317];
for ic = 1:length(cmr2)
    Planet.Cmeasured = cmr2(ic);
     for iS = 1:length(xS)
        Planet.xS = xS(iS);
        Seismic.mantleEOS = silicates{iS};
        Seismic.coreEOS = cores{iS};
       if rhos(ic,iS)
            disp([ic iS])
            Planet.rho_sil_withcore_kgm3 = rhos(ic,iS);
            PlanetProfile(Planet,Seismic,Params);
        end
    end
end

Planet.Ocean.comp='NaCl';
rmfield(Planet.Ocean,'fnTfreeze_K');
Params.LineStyle='-';
Params.CALC_NEW =1;
Params.CALC_NEW_REFPROFILES=1;
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY=0;
Params.HOLD = 0; % overlay previous run
Planet.FeCore=true;
Planet.Ocean.w_ocean_pct=2;  Planet.Tb_K = [266.38]; % pure water, temperatures at the bottom of the Ice Ih
% Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab'; % this needs to be substituted with the version that conserves Planet.rho_sil_withcore_kgm3 = 3518;
Seismic.mantleEOS = 'CI_anhydrous_differentiated_Ganymede_including_fluid_properties.tab';Planet.rho_sil_withcore_kgm3 = 3370;
% PlanetProfile(Planet,Seismic,Params)
rhos = [3350 3353 3351 3359 3361
        3328 3331 3330 3336 3338
        3305 3308 3305 3311 3313];
for ic = 1:length(cmr2)
    Planet.Cmeasured = cmr2(ic);
    for ifs = 1:length(xfes)
        Planet.xFeS = xfes(ifs);
        if rhos(ic,ifs)
            disp([ic ifs])
            Planet.rho_sil_withcore_kgm3 = rhos(ic,ifs);
            PlanetProfile(Planet,Seismic,Params);
        end
    end
end


