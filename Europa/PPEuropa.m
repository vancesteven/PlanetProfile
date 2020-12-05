function PPEuropa
%PPEuropa
Planet.name='Europa';
Planet.rho_kgm3 = 2989; % ±46 (Schubert et al. 2004)
Planet.R_m = 1561.0e3;% ±8.0 km
Planet.M_kg =4.7991e22;
Planet.gsurf_ms2 = 1.428; 
Planet.Tsurf_K = 110; 
Planet.Psurf_MPa = 0; 
Planet.Cmeasured = 0.346;
Planet.Cuncertainty = 0.005;
Planet.FeCore=true;
Planet.rhoFe = 8000; %8000
Planet.rhoFeS = 5150; %5150

%% salinities and temperatures at the bottom of the Ice Ih
% the vector of Tb needs to be monotonically increasing for the calculation
% of fluid electrical conductivities.

%%  Interior constraints imposed in Vance et al. 2014
% this information is not used here but set to keep PlanetProfile from
% throwing errors
mSi = 28.0855; mS = 32.065; mFe = 55.845; mMg = 24.305;
xOl = 0.44; % percentage of olivine - Javoy (1995) - Fe/Si = 0.909 Mg/Si = 0.531, Mg# = 0.369
xSi = (xOl+2*(1-xOl))*mSi/(xOl*184.295+(1-xOl)*244.3805); % mass fraction of sulfur in silicates
M_Earth_kg = 5.97e24;
xSiEarth = 0.1923; % Javoy in kg/kg in Exoplanets paper20052006-xSiSolar only in mole
xK = 1; %enrichment in K

%% Mantle Heat
%cold case  
Planet.kr_mantle = 4; % rock conductivity (Cammarano et al. 2006, Table 4)
Planet.Qmantle_Wm2 = 1.8e11/4/pi/Planet.R_m^2; % % this keeps the mantle temperature within the range that prevents melting.
% Planet.Qmantle_Wm2 = 2.2e11/4/pi/Planet.R_m^2; % this is more reasonable for radiogenic only
% Planet.Qmantle_Wm2 = 2.2e12/4/pi/Planet.R_m^2; % 
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
Seismic.QScore = 1e4;

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
%% set up the deeper interior properties. These are more petrologically consistent than simulations from Vance et al. 2018
Seismic.mantleEOS = 'Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.tab';
Seismic.coreEOS = '95Fe5S_core.tab';
xS = [0.05];
cmr2 = 0.346+[-0.005 0 0.005]; % Anderson et al. 1998
Planet.Cmeasured = cmr2(2);
Planet.xS = xS(1);
Planet.rho_sil_withcore_kgm3 = 3366;

% this information is not used here but set to keep PlanetProfile from
% throwing errors
Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
Planet.xFeS = 0.55; %0.25, mass fraction of sulfur in the core
Planet.xFe_core = 0.0279 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.0035; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.phi_surface = 0;


%% mgso4
Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
Params.LineStyle='--';
Params.wrefLine = '--';
Params.wref=[0 5 10 15];

Params.CALC_NEW =0;
Params.CALC_NEW_REFPROFILES=0;
Params.CALC_NEW_SOUNDSPEEDS=0;
Params.colororder = 'bm';

%10 wt%
Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [269.8 272.7]; 
PlanetProfile(Planet,Seismic,Params);

%1 wt%
Params.HOLD = 1; % overlay previous runs
Params.CALC_NEW =0;
Planet.Ocean.w_ocean_pct=1; Planet.Tb_K = [270.4  273.1]; % 265
PlanetProfile(Planet,Seismic,Params)

%% seawater
Params.HOLD = 1; % overlay previous runs
Planet.Ocean.comp='Seawater';
Params.LineStyle='-.';
Params.wref=[0 34];
Params.wrefLine = '-.';
Params.colororder = 'cb';

Params.CALC_NEW =0; % Set CALC_NEW options to 0 to re-use profile data when possible. It is recommended to keep CALC_NEW=1 except when intermediate parameters such as layer thicknesses will not change between runs.
Params.CALC_NEW_REFPROFILES=0;
Params.CALC_NEW_SOUNDSPEEDS=0;
Planet.Ocean.w_ocean_pct=gsw_SSO; Planet.Tb_K = [268.2 270.8 ];
PlanetProfile(Planet,Seismic,Params)

Params.CALC_NEW = 0; % Set CALC_NEW options to 0 to re-use profile data when possible. It is recommended to keep CALC_NEW=1 except when intermediate parameters such as layer thicknesses will not change between runs.
Planet.Ocean.w_ocean_pct=0.1*gsw_SSO; Planet.Tb_K = [270.0 272.5 ];
PlanetProfile(Planet,Seismic,Params)
