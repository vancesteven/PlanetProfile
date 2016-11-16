%PPTitan
 % this is the master program, and should be run from its containing
 % directory
Planet.name='Titan';
Planet.rho_kgm3 = 1879.8; %+/- 0.2, Jacobson et al. 2006
Planet.R_m = 2575e3; %+/- 2 km
Planet.M_kg =1.3452e23; % +/- 0.0002 
Planet.gsurf_ms2 = 1.428; 
Planet.Tsurf_K = 94; 
Planet.Psurf_MPa = 0; 
Planet.Cmeasured = 0.3438; % Fortes 2012, Iess 2010, 2012
% Planet.Cmeasured = 0.9*0.3438; % as suggested by Gao and Stevenson 2013
Planet.Cuncertainty = 0.0005;% Anderson et al. 2001 and Schubert et al. 2004 
Planet.FeCore=false;
Planet.rho_sil_withcore_kgm3 = 3400; %3250
Planet.xFeS = 0.25; %0.25
Planet.rhoFe = 8000; %8000
Planet.rhoFeS = 5150; %5150

Planet.Ocean.comp='MgSO4';
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [255 260 265 270]; % pure water, temperatures at the bottom of the Ice Ih
Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [252 255 260 265]; % 10 Wt% temperatures at the bottom of the Ice Ih

load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

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
Planet.kr_mantle = 4; % rock conductivity (Cammarano et al. 2006, Table 4)
Planet.Qmantle = 4e11; 
Planet.QHmantle = 0;

%% Seismic
Seismic.LOW_ICE_Q = 1; % divide Ice Q value by this number
Seismic.mantleEOS = 'ChondriteLL_Stx11.ext';
Seismic.QScore = 1e4;

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
Params.CALC_NEW =1;
Params.CALC_NEW_REFPROFILES=0;
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY=1;
Params.Pseafloor_MPa = 1000;
Params.nsteps_iceI = 100;
Params.nsteps_ocean = 450; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 100;
Params.nsteps_core = 10;
Params.savefigformat = 'epsc';
Params.wref=[0 5 10 15];
Params.colororder = 'mcbkgrm';
Params.Temps = [245 250 252.5 255 260 265 270 273];
colororder = Params.colororder(find(Params.Temps==Planet.Tb_K(1)):end);

%% Run the Calculation!
PlanetProfile(Planet,Seismic,Params)