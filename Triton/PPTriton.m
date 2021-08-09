%PPTriton_asymmetry
clear
Planet.name='Triton';
Params.cfg = config;
if Params.cfg.HOLD; clrAllProfiles; clrAllLayered(Planet.name); end
%% &&& Orbital and plotting parameters for use in LayeredInductionResponse
Planet.ionos_bounds = [250e3 450e3];
Planet.ionosPedersen_sig = [1e-16 0.05];

Planet.R_m = 1353.4e3; %±0.9km "Planetary Satellite Physical Parameters". JPL (Solar System Dynamics). Archived from the original on January 18, 2010. Retrieved October 26, 2011.
Planet.M_kg =2.14e22; %± "Planetary Satellite Physical Parameters". JPL (Solar System Dynamics). Archived from the original on January 18, 2010. Retrieved October 26, 2011.
Planet.Tsurf_K = 38; 
Planet.Psurf_MPa = 0;
Planet.Cmeasured = 0.315; %  Hussmann 2006 estimate for *Pluto* (not actually measured)
Planet.Cuncertainty = 0.008;% 
Planet.FeCore=false;
    Planet.xFe_core = 0.75;
    Planet.xFeS = 0.25; %0.25
    Planet.rhoFe = 8000; %8000
    Planet.rhoFeS = 5150; %5150
Planet.rho_sil_withcore_kgm3 = 2400; % Iess et al. 2014

% WARNING: The following line was copied from PPCallisto.m because it is
% required by PlanetProfile.m in the current version and does not appear in
% this file. An issue has been opened on GitHub. Delete this comment when
% the value has been corrected.
Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model

%%  Interior constraints imposed in Vance et al. 2014
mSi = 28.0855; mS = 32.065; mFe = 55.845; mMg = 24.305;
xOl = 0.44; % percentage of olivine - Javoy (1995) - Fe/Si = 0.909 Mg/Si = 0.531, Mg# = 0.369
xSi = (xOl+2*(1-xOl))*mSi/(xOl*184.295+(1-xOl)*244.3805); % mass fraction of sulfur in silicates
M_Earth_kg = 5.97e24;
xSiEarth = 0.1923; % Javoy in kg/kg in Exoplanets paper20052006-xSiSolar only in mole
xK = 1; %enrichment in K
Hrad0 = 24e12*xSi/xSiEarth/M_Earth_kg;

%% Mantle Heat
%cold case  
Planet.kr_mantle = 4; % rock conductivity (Cammarano et al. 2006, Table 4)
Planet.EQUIL_Q = 0;
Planet.Qmantle_Wm2 = 2.7e8/4/pi/Planet.R_m^2; % Chen et al. 2014
Planet.QHmantle = 0;

%% Seismic
Seismic.LOW_ICE_Q = 1; % divide Ice Q value by this number

Planet.POROUS_ROCK = 1;
Planet.phi_surface = 0.8;

Seismic.QScore = 1e4;
Seismic.SMOOTH_VROCK = 5; % smooth over N neighboring rows and columns in vp and vs

%Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab';
Seismic.mantleEOS =  'CI_hhph_DEW17_nofluid_nomelt_685.tab';
%Seismic.mantleEOS =  'Comet_67P_CG_v7_excluding_fluid_properties.tab';
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
Params.savefigformat = 'epsc';
Params.foursubplots =1;
Params.Legend = false;
Params.LegendPosition = 'North';
Params.ylim = [910 1170];
Params.Pseafloor_MPa = 500; % be careful to make sure this is in a reasonable range
Params.nsteps_iceI = 200;
Params.nsteps_ocean = 100; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 1500;
Params.nsteps_core = 100;
Params.Temps = [245 250 252.5 255 260 265 270];

%% Run the Calculation!
% 
Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

Params.LineStyle='--';
Params.wref=[0 5 10 15];
Params.wrefLine = '--';
Params.colororder = 'cmb';

Planet.Ocean.w_ocean_pct=1; Planet.Tb_K = [266.0]; % 112 km thick ice

outPlanet = PlanetProfile(Planet,Seismic,Params);
