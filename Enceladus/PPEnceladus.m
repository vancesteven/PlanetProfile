%PPEnceladus
clear
Planet.name='Enceladus';

Params.cfg = config;
if Params.cfg.HOLD; clrAllProfiles; end

%% &&& Bulk and surface properties
Planet.rho_kgm3 = 1610; % Thomas 2010
Planet.R_m = 252.1e3; %± 200 Thomas  2010
Planet.M_kg =1.08022e20; 
Planet.gsurf_ms2 = 0.113; 
Planet.Tsurf_K = 75; 
Planet.Psurf_MPa = 0; 
Planet.Cmeasured = 0.335; % Iess et al. 2014
Planet.Cuncertainty = 0.001;% 
Planet.FeCore=false;
    Planet.xFeS = 0.25; %0.25
    Planet.rhoFe = 8000; %8000
    Planet.rhoFeS = 5150; %5150
% Planet.rho_sil_withcore_kgm3 = 2400; % Iess et al. 2014
Planet.rho_sil_withcore_kgm3 = 2700; 

% WARNING: The following line was copied from PPCallisto.m because it is
% required by PlanetProfile.m in the current version and does not appear in
% this file. An issue has been opened on GitHub. Delete this comment when
% the value has been corrected.
Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.POROUS_ICE = 0;
Planet.phi_surface = 0.5;

% Planet.Ocean.comp='NH3';
% % load L_IceNH3.mat
% Planet.Ocean.w_ocean_pct=5;
% Planet.Tb_K = [245 260 265]; % 5 Wt% temperatures at the bottom of the Ice Ih Van dishoek, cited by Engel et al. 1994 for Titan
% Tb = [251.5 260 270]; % 1 Wt% temperatures at the bottom of the Ice Ih
% Van dishoek, cited by Engel et al. 1994 for Titan


%%  Interior constraints imposed in Vance et al. 2014
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
Planet.kr_mantle = 4; % rock conductivity (Cammarano et al. 2006, Table 4)
Planet.EQUIL_Q = 0;
% Planet.Qmantle_Wm2 = 2.7e8/4/pi/Planet.R_m^2; % Chen et al. 2014
 Planet.Qmantle_Wm2 = 2.7e4/4/pi/Planet.R_m^2; % kluge
% Planet.Qmantle_Wm2 = 16e9/4/pi/Planet.R_m^2; % Howett et al. 2011
Planet.QHmantle = 0;
%hot case Qm = 2.1e11+8.5e11; %W

%% Seismic
Seismic.LOW_ICE_Q = 1; % divide Ice Q value by this number
% Seismic.mantleEOS = 'ChondriteLL_Stx11.ext'; % this was the original
% chondrite function built by svance

% this is the porous model described in the paper 2017:
Planet.POROUS_ROCK = 1;
Seismic.mantleEOS = 'echon_678_1.tab'; Planet.phi_surface = 0.8;% 
% Seismic.mantleEOS = 'pyrohp_sat_678_1.tab'; Planet.phi_surface = 0.5;% this uses the procedure implemented by F. Cammarano

% % this is the hydrated model described in the paper 2017:
% Planet.POROUS_ROCK = 0;
% Seismic.mantleEOS = 'pyrohp_sat_678_1.tab'; % this uses the procedure implemented by F. Cammarano
Seismic.mantleEOS = 'echon_hp_sat_PX678_14GPa.tab'; % this uses the procedure implemented by F. Cammarano


% Seismic.mantleEOS = 'CM_hhph_DEW17_fluid_nomelt_685.tab';
Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab';

% Seismic.mantleEOS = 'echonhp_sat_1.tab'; % this uses the procedure implemented by F. Cammarano

Seismic.QScore = 1e4;
Seismic.SMOOTH_VROCK = 5; % smooth over N neighboring rows and columns in vp and vs


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
Params.Legend = false;
Params.LegendPosition = 'North';
Params.ylim = [910 1170];
Params.Pseafloor_MPa = 20;
Params.nsteps_iceI = 20;
Params.nsteps_ocean = 100; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 1500;
Params.nsteps_core = 100;
Params.Temps = [245 250 252.5 255 260 265 270];


%% Run the Calculation!
% Planet.Ocean.w_ocean_pct=10;Planet.Tb_K = [272.8 272.9 273 273.1]; % pure water, temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
% 
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [273.1 273.15]; % pure water, temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
% 
Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
Params.LineStyle='--';
Params.colororder = 'cm';
Params.wrefLine = '--';
Params.wref=[0 5 10 15];
% 
Planet.Ocean.w_ocean_pct=5; Planet.Tb_K = [272.72 273.12];
%outPlanet = PlanetProfile(Planet,Seismic,Params);

% running pure water for the MgSO4 case illustrates >1oC error in the Margules parameterization
Params.LineStyle='-';
Params.colororder = 'cm';
Planet.Ocean.w_ocean_pct=gsw_SSO/3;  Planet.Tb_K = 271.0; % pure water, 
outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
% Planet.Ocean.comp='Seawater';
% Params.LineStyle='-.';
% Params.wref=[0 34 68];
% Params.wrefLine = '-.';
% Params.colororder = 'cm';
% 
% Planet.Ocean.w_ocean_pct=gsw_SSO; Planet.Tb_K = [270.82  271.16];
% outPlanet = PlanetProfile(Planet,Seismic,Params);

% Params.LineStyle='-';
% Params.colororder = 'cm';
% Planet.Ocean.w_ocean_pct=0;  Planet.Tb_K =  [272.74 273.08]; % pure water, 
% outPlanet = PlanetProfile(Planet,Seismic,Params);

Params.wrefLine =  ':';
Params.wref=[3 5 10];
Planet.Ocean.comp='NH3';
load L_IceNH3_DATA.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

Params.LineStyle =  ':';
Planet.Ocean.w_ocean_pct=3; Planet.Tb_K = [269.535 269.905]; % 0 Wt% temperatures at the bottom of the Ice Ih
%outPlanet = PlanetProfile(Planet,Seismic,Params);