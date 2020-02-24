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



%Planet.Tb_K = [250 255 260 265 270]; % 15 Wt% temperatures at the bottom of the Ice Ih
% %Planet.Tb_K = [252.5 255 260 265 270]; %3 and 5 Wt% temperatures at the bottom of the Ice Ih
%  Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [252 255 260 265 270]; % 10 Wt% temperatures at the bottom of the Ice Ih
Params.LineStyle='--';
% Seismic.mantleEOS = 'chon_678_1.tab'; %(3500)
% Seismic.mantleEOS = 'epyro_678_1.tab'; %(3450)
% Seismic.mantleEOS = 'pyrohp_sat_678_1.tab'; %(3305)
% Seismic.mantleEOS = 'epyro_678_1.tab'; % (3450) this uses the procedure implemented by F. Cammarano
% Seismic.mantleEOS = 'chonhp_sat_678.tab';% (3450)

Seismic.mantleEOS = 'echon_hp_sat_PX678_14GPa.tab';% (3450)
Seismic.mantleEOS = 'pyro_678_1.tab'; % (3430) this uses the procedure implemented by F. Cammarano

% Seismic.mantleEOS = 'CV3hy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.55; %0.25
% Planet.xFe_core = 0.0279 ; % this is the total Fe  in Fe and FeS
% Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3740; 

% Seismic.mantleEOS = 'CIhy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Planet.xFeS_meteoritic = 0.0908; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.2; %0.25
% Planet.xFe_core = 0.0583 ; % this is the total Fe  in Fe and FeS
% Planet.XH2O = 0.169; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3440;


Seismic.mantleEOS = 'CM2hy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
Planet.xFeS_meteoritic = 0.0676; %CM2 mean from Jarosewich 1990
Planet.xFeS = 0.2; %0.25
Planet.xFe_core = 0.0463 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.rho_sil_withcore_kgm3 = 3730;

interiors = {
% 'CI_hhph_DEW17_fluid_nomelt_685.tab'  
'CI_hhph_DEW17_nofluid_nomelt_685.tab'
% 'CM_hhph_DEW17_fluid_nomelt_685.tab'  
'CM_hhph_DEW17_nofluid_nomelt_685.tab'
% 'CV_hhph_DEW17_fluid_nomelt_685.tab'  
'CV_hhph_DEW17_nofluid_nomelt_685.tab'
% 'H_hhph_DEW17_685_fluid_nomelt.tab'   
'H_hhph_DEW17_685_nofluid_nomelt.tab'   
% 'L_hhph_DEW17_685_v1_fluid_nomelt.tab'  
'L_hhph_DEW17_685_v1_nofluid_nomelt.tab'
'Simple_CI_HS_green_PP.tab'
'Simple_CM_HS_green_PP.tab'
'Simple_CV_HS_green_PP.tab'
};


% [125 100 75 50]
Params.NOPLOTS = 0; %allows user to limit recreating plots & figures after each run

% Planet.Cmeasured = 0.3115;
cmr2 = 0.3115+[-0.0028 0 0.0028];
xfes = [0 .05 .01 0.15 0.20];
Planet.Cmeasured = 0.3115;
Planet.xFeS = 0.2; %0.25

Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
% rmfield(Planet.Ocean,'fnTfreeze_K');
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

Params.CALC_NEW =1; % Set CALC_NEW options to 0 to re-use profile data when possible. It is recommended to keep CALC_NEW=1 except when intermediate parameters such as layer thicknesses will not change between runs.
Params.CALC_NEW_REFPROFILES=1;
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY=0;
Planet.FeCore=true;
% Seismic.coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab'; 

% Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab'; Planet.rho_sil_withcore_kgm3 = 3505; % this needs to be substituted with
% the version that partitions Fe and S between the core and silicates
% Seismic.mantleEOS = 'CI_hydrous_differentiated_Ganymede_including_fluid_properties.tab';
Seismic.mantleEOS = 'CI_anhydrous_differentiated_Ganymede_including_fluid_properties.tab';
% Seismic.mantleEOS = 'CI_hydrous_differentiated_Ganymede_excluding_fluid_properties.tab';
Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [254.75 259.05 263.1 266]; % 10 Wt% temperatures at the bottom of the Ice Ih
% rhos = [3378 3378 3381 3571 3388
%         3357 3360 3357 3367 3367
%         3335 3335 3337 3342 3344];
% for ic = 1:length(cmr2)
%     Planet.Cmeasured = cmr2(ic);
%     for ifs = 1:length(xfes)
%         Planet.xFeS = xfes(ifs);
%         if rhos(ic,ifs)
%             disp([ic ifs])
%             Planet.rho_sil_withcore_kgm3 = rhos(ic,ifs);
%             PlanetProfile(Planet,Seismic,Params);
%         end
%     end
% end

Seismic.mantleEOS = 'CM_anhydrous_differentiated_Ganymede_including_fluid_properties.tab';
% Seismic.mantleEOS = 'CM_hydrous_differentiated_Ganymede_including_fluid_properties.tab';
% Seismic.mantleEOS = 'CM_hydrous_differentiated_Ganymede_excluding_fluid_properties.tab';
rhos = [3378 3378 3381 3571 3388
        3357 3360 3357 3367 3367
        3335 3335 3337 3342 3344];
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



Params.CALC_NEW =1;
Params.CALC_NEW_REFPROFILES=1;
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY=0;
Params.HOLD = 0; % overlay previous run
Planet.FeCore=true;
% Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab'; % this needs to be substituted with the version that conserves Planet.rho_sil_withcore_kgm3 = 3518;
Seismic.mantleEOS = 'CI_anhydrous_differentiated_Ganymede_including_fluid_properties.tab';Planet.rho_sil_withcore_kgm3 = 3370;
Planet.Ocean.w_ocean_pct=0;  Planet.Tb_K = [256.5  260.6  264.4 267.7]; % pure water, temperatures at the bottom of the Ice Ih
Params.LineStyle='-';
% PlanetProfile(Planet,Seismic,Params)
rhos = [3348 3351 3349 3357 3359
        3327 3329 3327 3334 3335
        3304 3307 3304 3309 3313];
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


Params.CALC_NEW =1;
Params.CALC_NEW_REFPROFILES=1;
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
    for ifs = 1:length(xfes)
        Planet.xFeS = xfes(ifs);
        if rhos(ic,ifs)
            disp([ic ifs])
            Planet.rho_sil_withcore_kgm3 = rhos(ic,ifs);
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


