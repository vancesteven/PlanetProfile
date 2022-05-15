%PPTitan
clear
Planet.name='Titan';
  
Params.cfg = config;
if Params.cfg.HOLD; clrAllProfiles; end

%% &&& Bulk and surface properties
% Planet.rho_kgm3 = 1879.8; %+/- 0.2, Jacobson et al. 2006
Planet.R_m= 2574.7e3; % christophe's value
Planet.M_kg = 1.3453e23; %christophe's value
% Planet.R_m = 2575e3; %+/- 2 km
% Planet.M_kg =1.3452e23; % +/- 0.0002 
Planet.Tsurf_K = 94; 
Planet.Psurf_MPa = 0; 
Planet.xFeS = 0; %0.25
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
Planet.Qmantle_Wm2 = 1.2e10/4/pi/Planet.R_m^2; 
% Planet.Qmantle_Wm2 = 4e12/4/pi/Planet.R_m^2; 
Planet.QHmantle = 0;

%% Porosity of the rock
Planet.POROUS_ROCK = 0; % porosity makes no difference for Titan because the seafloor pressures are too high
Planet.phi_surface = 1;

%% Seismic
Seismic.LOW_ICE_Q = 1; % divide Ice Q value by this number
% Seismic.mantleEOS = 'echon_hp_sat_PX678_14GPa.tab'; % this uses the procedure implemented by F. Cammarano; this includes Ks and Gs. I had to rerun perlex (6.6.3). not sure why
% Seismic.mantleEOS = 'chonhp_sat_678.tab';% (3300)

% Seismic.mantleEOS = 'pyrohp_sat_678_1.tab'; %  (3000) this uses the procedure implemented by F. Cammarano; this includes Ks and Gs. I had to rerun perlex (6.6.3). not sure why
% Seismic.mantleEOSname = 'pyrohpsat';

% Seismic.mantleEOS = 'CV_hhph_DEW17_nofluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CV_hhphD17nofluid';

Seismic.mantleEOS = 'CV_hhph_DEW17_fluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CV_hhphD17fluid';

% Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CM_hhphD17nofluid';

% Seismic.mantleEOS = 'CM_hhph_DEW17_fluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CM_hhphD17fluid';

% Seismic.mantleEOS = 'CI_hhph_DEW17_nofluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CI_hhphD17nofluid';

% Seismic.mantleEOS = 'CI_hhph_DEW17_fluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CI_hhphD17fluid';

Seismic.QScore = 1e4;
Seismic.SMOOTH_VROCK = 1; % smooth over N neighboring rows and columns in vp and vs

%Attenuation Parameters Based on those Described in Cammarano et al. 2006
% ice I
Seismic.B_aniso_iceI = 0.56;
Seismic.gamma_aniso_iceI = 0.2;
Seismic.g_aniso_iceI = 22; %C2006
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
Params.Pseafloor_MPa = 1800;
Params.Legend = false;
Params.LegendPosition = 'southeast';
Params.ylim = [925 1350];
Params.nsteps_iceI = 100;
Params.nsteps_ocean = 450; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 100;
Params.nsteps_core = 10;

Params.colororder = 'cbmkgrm';
Params.Temps = [245 250 252.5 255 260 265 270 273];

%% Run the Calculation!
Planet.Cuncertainty = 0.0005;%
Planet.Cmeasured = 0.3438; % Fortes 2012, Iess 2010, 2012
% Planet.Cmeasured = 0.3318; % in prep Sotin

Planet.FeCore=false;
% 
% Planet.FeCore=true;
% Planet.rho_sil_withcore_kgm3 = 2500;

% placeholder so Matlab doesn't complain about missing fields
Planet.xFeS_meteoritic = 0.0676; %CM2 mean from Jarosewich 1990
Planet.xFeS = 1; %0.25
Planet.xFe_core = 0.0463 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model

% Comparison of MgSO4 EOS pure water values is close to values for ammonia
% EOS.  There are diffences in both the melting temperatures and fluid
% thermodynamics

Params.LineStyle =  '--';
Params.wrefLine =  '--';
Params.wref=[0 5 10 15];
Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [255 260 265 270]; % 0 Wt% temperatures at the bottom of the Ice Ih
% Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [252 255 260 266]; % 10 Wt% temperatures at the bottom of the Ice Ih
Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = 266; % 10 Wt% temperatures at the bottom of the Ice Ih; as currently in the manuscript
% Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [251.1 262 266]; % 10 Wt% temperatures at the bottom of the Ice Ih
outPlanet = PlanetProfile(Planet,Seismic,Params);

% %== Supporting pure NaCl oceans is in development.
% Planet.Ocean.comp='NaCl';
% Params.LineStyle =  ':';
% Params.wrefLine =  ':';
% Params.wref=[0 5 10 15];
% load L_Ice_MgSO4.mat
% Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [255 260 265 270]; % 0 Wt% temperatures at the bottom of the Ice Ih
% Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [252 255 260 266]; % 10 Wt% temperatures at the bottom of the Ice Ih
%Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [246.9 260]; % 10 Wt% temperatures at the bottom of the Ice Ih; as currently in the manuscript
% Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [251.1 262 266]; % 10 Wt% temperatures at the bottom of the Ice Ih
%outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
%==
% 
% Params.wrefLine =  '-.';
% Params.wref=[3 5 10];
% Planet.Ocean.comp='NH3';
% load L_IceNH3_DATA.mat
% Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
% 
% Params.LineStyle =  '-';
% % Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [255 260 265 270]; % 0 Wt% temperatures at the bottom of the Ice Ih
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [255 265 268]; % 0 Wt% temperatures at the bottom of the Ice Ih; for the paper
% % Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [254 265 268]; % 0 Wt% temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
% Params.LineStyle =  '-.';
% 
% Planet.Ocean.w_ocean_pct=3; Planet.Tb_K = [250 260 264]; % 3 Wt% temperatures at the bottom of the Ice Ih; as currently in the manuscript, thicknest ice is 150 km
% % Planet.Ocean.w_ocean_pct=3; Planet.Tb_K = [249 260 264]; % 3 Wt% temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
% Planet.Ocean.w_ocean_pct=7; Planet.Tb_K = [250 260]; % 7 Wt% temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);

% Params.wrefine =  '--';
% Params.LineStyle =  '-';
% 
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [255 260 265 270]; % pure water, temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% 

% TESTING THE INFLUENCE OF A LOWER MOMENT OF INERTIA
% 
% Planet.Cmeasured = 0.31; % 0.9*0.3438;  as suggested by Gao and Stevenson 2013
% Planet.FeCore=true;
% Planet.rho_sil_withcore_kgm3 = 3400; %3250
% 
% Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [252 255 260 265]; % 10 Wt% temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [255 260 265 270]; % pure water, temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);
