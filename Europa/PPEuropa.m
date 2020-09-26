function PPEuropa
%PPEuropa
Planet.name='Europa';
Planet.rho_kgm3 = 2989; % �46 (Schubert et al. 2004)
%note: Schubert et al. 2004 cite the Anderson C/MR2 as 0.3115�0.0028.  This
%is incorrect, as the value cited everywhere else is consistent with the
%Anderson et al. (1996) value of C/MR2=0.3105\pm0.0028 used here
%MMD note July 17 2018: But Anderson et al. 1998 (Science) reported a
%preferred value of C/MR2=0.346, which is actually used here.
Planet.R_m = 1561.0e3;% �8.0 km
Planet.M_kg =4.7991e22;
Planet.gsurf_ms2 = 1.428; 
Planet.Tsurf_K = 110; 
Planet.Psurf_MPa = 0; 
Planet.Cmeasured = 0.346;
Planet.Cuncertainty = 0.005;
%note: Schubert et al. 2004 cite the Anderson C/MR2 as 0.3115�0.0028.  This
%is incorrect, as the value cited everywhere else is consistent with the
%Anderson et al. (1996) value of C/MR2=0.3105\pm0.0028 used here
%MMD note July 17 2018: But Anderson et al. 1998 (Science) reported a
%preferred value of C/MR2=0.346, which is actually used here.
Planet.FeCore=true;
Planet.rhoFe = 8000; %8000
Planet.rhoFeS = 5150; %5150
%Planet.rhoPoFeFCC = 5455; %�40; WIP July 17 2018; Density of pyrrhottite plus 
%face-centered cubic iron predicted to stable allowing the maximum amount 
%of sulfur in Europa's current core without crossing the solidus, according 
%to our calculation using the Saxena and Eriksson 2015 (CALPHAD) EoS, 
%modified by Eleanor Green and Jamie Connolly for PerpleX 6.8.3.

%% salinities and temperatures at the bottom of the Ice Ih
% the vector of Tb needs to be monotonically increasing for the calculation
% of fluid electrical conductivities.

%%  Interior constraints imposed in Vance et al. 2014
mSi = 28.0855; mS = 32.065; mFe = 55.845; mMg = 24.305;
xOl = 0.44; % percentage of olivine - Javoy (1995) - Fe/Si = 0.909 Mg/Si = 0.531, Mg# = 0.369
%mOl = 2*((1-0.369)*58.85+0.369*24.31)+28.0855+4*16=184.295
%mPx = 2*((1-0.369)*58.85+0.369*24.31+28.0855+3*16) =244.3805
xSi = (xOl+2*(1-xOl))*mSi/(xOl*184.295+(1-xOl)*244.3805); % mass fraction of sulfur in silicates
M_Earth_kg = 5.97e24;
xSiEarth = 0.1923; % Javoy in kg/kg in Exoplanets paper20052006-xSiSolar only in mole
xK = 1; %enrichment in K

%% Mantle Heat
%cold case  
Planet.kr_mantle = 4; % rock conductivity (Cammarano et al. 2006, Table 4)
% Planet.Qmantle_Wm2 = 1e11/4/pi/Planet.R_m^2; %
Planet.Qmantle_Wm2 = 2.2e11/4/pi/Planet.R_m^2; % this is more reasonable for radiogenic only
Planet.QHmantle = 0;
Planet.EQUIL_Q = 0;
%hot case Qm = 2.1e11+8.5e11; %W
% Qmantle = 1.3e11; 
% QHmantle = 8.5e11;

%% Porosity of the rock
Planet.POROUS_ROCK = 0;
% Planet.PEFF =0;

%% Seismic
Seismic.LOW_ICE_Q = 1; % divide Ice Q value by this number
% Seismic.mantleEOS = 'epyrohp_sat_1.tab'; % low density (2700) this uses the procedure implemented by F. Cammarano


% Seismic.mantleEOS = 'chon_678_1.tab'; %(3500)
% Seismic.mantleEOS = 'pyrohp_sat_678_1.tab'; %(2820)
% Seismic.mantleEOS = 'epyro_678_1.tab'; %(3450)
% Seismic.mantleEOS = 'echon_hp_sat_PX678_14GPa.tab';% (3100)
% Seismic.mantleEOS = 'pyrohy_678v2_1.tab'; %(3420)
% Seismic.mantleEOS = 'echonhy1wt_678_1.tab'; %(3459)


% Seismic.mantleEOS = 'pyrohy_1wt_14GPa.tab'; % (3300) this uses the procedure implemented by F. Cammarano
% Seismic.mantleEOS = 'echonhp_1wt_6GPa.tab'; % too dense (3400) this uses the procedure implemented by F. Cammarano
%  Seismic.mantleEOS = 'epyro_1.tab'; % too dense (3400) this uses the procedure implemented by F. Cammarano
% Seismic.mantleEOS = 'echon_1.tab'; % too dense (3500) this uses the procedure implemented by F. Cammarano
% Seismic.mantleEOS = 'chon_1.tab'; % too dense (3500) this uses the procedure implemented by F. Cammarano
Seismic.QScore = 1e4;

Seismic.coreEOS = 'sulfur_core_partition_SE15_1pctSulfur.tab';

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

%% Global Model Parameters
Params.HOLD = 0; % overlay previous runs
Params.foursubplots =1;
Params.Legend = 0;
Params.LegendPosition = 'North'; 
Params.ylim = [910 1230];
Params.Pseafloor_MPa = 350;
Params.nsteps_iceI = 200;
Params.nsteps_ocean = 350; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 500;
Params.nsteps_core = 10;
Params.savefigformat = 'epsc';
Params.colororder = 'mcbkgrm';
Params.Temps = [250 252.5 255 260 265 270 273];
Params.NOPLOTS = 0; %allows user to control recreating & display of plots/figures after each run

%% Run the Calculation!
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY = 1;
% 
Seismic.mantleEOS = 'chon_678_1.tab'; %(3440) % this did not exclude nasGL and faGL and so had many nan entries
Planet.rho_sil_withcore_kgm3 = 3539;
% Seismic.mantleEOS = 'pyrohy_678v2_1.tab'; %(3440) % this did not exclude nasGL and faGL and so had many nan entries
% Planet.rho_sil_withcore_kgm3 = 3425;
% Seismic.mantleEOS = 'chonhp_sat_678.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Planet.rho_sil_withcore_kgm3 = 2975;
% Seismic.mantleEOS = 'Simple_CI_HS_green_PP.tab';% CI chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
% Planet.rho_sil_withcore_kgm3 = 2975;
% Seismic.mantleEOS = 'Simple_CM_HS_green_PP.tab';% CM chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
% Planet.rho_sil_withcore_kgm3 = 2975;
% Seismic.mantleEOS = 'Simple_CV_HS_green_PP.tab';% CV chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
% Planet.rho_sil_withcore_kgm3 = 2975;

Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
Params.LineStyle='--';
Params.wrefLine = '--';
Params.wref=[0 5 10 15];

Params.CALC_NEW =1;
Params.CALC_NEW_REFPROFILES=1;
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.colororder = 'cm';
Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [269.8  272.7]; % 265


Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
Planet.xFeS = 0.55; %0.25, mass fraction of sulfur in the core
Planet.xFe_core = 0.0279 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.rho_sil_withcore_kgm3 = 3644;
Planet.phi_surface = 0;

PlanetProfile(Planet,Seismic,Params)

Params.HOLD = 1; % overlay previous runs
% Params.LineStyle='-';
% Params.colororder = 'cm';
% Planet.Ocean.w_ocean_pct=0;  Planet.Tb_K =  [270.4 273.1]; % pure water,  265.7
% PlanetProfile(Planet,Seismic,Params)

Planet.Ocean.comp='Seawater';
Params.LineStyle='-.';
Params.wref=[0 34];
Params.wrefLine = '-.';
Params.colororder = 'cm';

Params.CALC_NEW =1; % Set CALC_NEW options to 0 to re-use profile data when possible. It is recommended to keep CALC_NEW=1 except when intermediate parameters such as layer thicknesses will not change between runs.
Params.CALC_NEW_REFPROFILES=1;
Params.CALC_NEW_SOUNDSPEEDS=1;

Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
Planet.xFeS = 0.55; %0.25, mass fraction of sulfur in the core
Planet.xFe_core = 0.0279 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.rho_sil_withcore_kgm3 = 3644;
Planet.phi_surface = 0;

Planet.Ocean.w_ocean_pct=gsw_SSO; Planet.Tb_K = [268.2 270.8 ];
 
% Seismic.mantleEOS = 'chon_678_1.tab'; %(3440) % this did not exclude nasGL and faGL and so had many nan entries
% Planet.rho_sil_withcore_kgm3 = 3539;
% Planet.phi_surface = 0;
% PlanetProfile(Planet,Seismic,Params)
% Params.HOLD = 1;
% % Planet.PEFF =1;
% % PlanetProfile(Planet,Seismic,Params)
% % Planet.PEFF =0;
% % 
% Planet.phi_surface = 0.8;
% PlanetProfile(Planet,Seismic,Params)
% % Planet.PEFF =1;
% % PlanetProfile(Planet,Seismic,Params)
% % Planet.PEFF =0;

% Seismic.mantleEOS = 'chonhp_sat_678.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Seismic.mantleEOS = 'CI_expanded_678.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% 
% Seismic.mantleEOS = 'CM2hy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Planet.xFeS_meteoritic = 0.0676; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.2; %0.25, mass fraction of sulfur in the core
% Planet.xFe_core = 0.0463 ; % this is the total Fe  in Fe and FeS
% Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3630;
% Planet.phi_surface = 0;
% PlanetProfile(Planet,Seismic,Params)

% Planet.PEFF =1;
% PlanetProfile(Planet,Seismic,Params)
% Planet.PEFF =0;
% 
% Planet.phi_surface = 0.8;
% PlanetProfile(Planet,Seismic,Params)
% % Planet.PEFF =1;
% % PlanetProfile(Planet,Seismic,Params)
% 
Seismic.mantleEOS = 'CV3hy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
Planet.xFeS = 0.55; %0.25, mass fraction of sulfur in the core
Planet.xFe_core = 0.0279 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.rho_sil_withcore_kgm3 = 3644;
Planet.phi_surface = 0;
PlanetProfile(Planet,Seismic,Params)
% 
% Seismic.mantleEOS = 'CIhy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Planet.xFeS_meteoritic = 0.0908; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.2; %0.25, mass fraction of sulfur in the core
% Planet.xFe_core = 0.0583 ; % this is the total Fe  in Fe and FeS
% Planet.XH2O = 0.169; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3644;
% Planet.phi_surface = 0;
% PlanetProfile(Planet,Seismic,Params)
% 
% Seismic.mantleEOS = 'Simple_CI_HS_green_PP.tab';% CI chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
% Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.0; %mass fraction of sulfur in the core
% Planet.xFe_core = 0.0279 ; % this is the total Fe in Fe and FeS
% Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3644;
% Planet.phi_surface = 0;
% PlanetProfile(Planet,Seismic,Params)

% Seismic.mantleEOS = 'Simple_CM_HS_green_PP.tab';% CM chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
% Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.0; %mass fraction of sulfur in the core
% Planet.xFe_core = 0.0279 ; % this is the total Fe in Fe and FeS
% Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3644;
% Planet.phi_surface = 0;
% PlanetProfile(Planet,Seismic,Params)

% Seismic.mantleEOS = 'Simple_CV_HS_green_PP.tab';% CV chondrite material minus Fe core, computed with Green et al. 2016 (JMG) solution model and Lodders and Fegley 1998
% Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.0; %mass fraction of sulfur in the core
% Planet.xFe_core = 0.0279 ; % this is the total Fe in Fe and FeS
% Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3644;
% Planet.phi_surface = 0;
% PlanetProfile(Planet,Seismic,Params)