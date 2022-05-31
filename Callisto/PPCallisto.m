%PPCallisto
clear
Planet.name='Callisto';

Params.cfg = config;
if Params.cfg.HOLD; clrAllProfiles; clrAllLayered(Planet.name); end
%% &&& Orbital and plotting parameters for use in LayeredInductionResponse
Planet.peaks_Hz = [5.4584e-05 2.7294e-5 6.892e-7];
Planet.f_orb = 2*pi/17/86400; % radians per second
Params.wlims = [log(0.01) log(100000)];
% Get Fourier spectrum data
try 
    load(['FTdata' Planet.name],'FTdata');
catch
    dlNow = FTdataMissingWarning();
    if dlNow; load(['FTdata' Planet.name],'FTdata'); end
end
Planet.ionos_bounds = 100e3;
Planet.ionosPedersen_sig = 30/100e3;
Planet.ionos_only = [];
Planet.PLOT_SIGS = true;
Planet.ADD_TRANSITION_BOUNDS = false;
 
%% &&& Bulk and surface properties
Planet.rho_kgm3 = 1834.4; % Schubert et al. 2004, Anderson et al. 2001; ±3.4
Planet.R_m = 2410.3e3; %±1.5e3
Planet.M_kg =1.0759e23; 
Planet.gsurf_ms2 = 1.428; 
Planet.Tsurf_K = 110; 
Planet.Psurf_MPa = 0; 
Planet.FeCore=false;
Planet.xFeS = 0; %0.25
Planet.rhoFe = 8000; %8000
Planet.rhoFeS = 5150; %5150

Planet.Ocean.comp='MgSO4';
%Planet.Tb_K = [252.5 255 260 265 270]; %3 and 5 Wt% temperatures at the bottom of the Ice Ih
%shell for the 3 Wt% case

load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

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
Planet.EQUIL_Q = 0;
Planet.kr_mantle = 4; % rock conductivity (Cammarano et al. 2006, Table 4)
% Planet.Qmantle_Wm2 = 1.1e8/4/pi/Planet.R_m^2; % this also works for the
% saturated pyrolite, and marginally for the saturated chondrite
Planet.Qmantle_Wm2 = 1.3e11/4/pi/Planet.R_m^2; % this works for the saturated pyrolite
% Planet.Qmantle_Wm2 = 1.1e12/4/pi/Planet.R_m^2; % test
Planet.QHmantle = 0;

%% Porosity of the rock
Planet.POROUS_ROCK = 0;


%% Seismic
Seismic.LOW_ICE_Q = 1; % divide Ice Q value by this number
Seismic.SMOOTH_VROCK = 1; % smooth over N neighboring rows and columns in vp and vs
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
Params.foursubplots =1;
Params.Legend =0;
Params.LegendPosition = 'north';
Params.ylim = [925 1350];
Params.Pseafloor_MPa = 1000;
Params.nsteps_iceI = 100;
Params.nsteps_ocean = 450; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 100;
Params.nsteps_core = 10;
Params.wref=[0 5 10 15];
Params.colororder = 'cbmkgrm';
Params.Temps = [245 250 252.5 255 260 265 270 273];


%% Run the Calculation!
Planet.Cmeasured = 0.3549; 
Planet.Cuncertainty = 0.0042;% Anderson et al. 2001 and Schubert et al. 2004 
Seismic.mantleEOS = 'pyrohp_sat_678_1.tab';% this was used in the published paper (Vance et al. JGR 2017)
% Seismic.mantleEOS = 'chonhp_sat_678.tab';% (3300)

Planet.FeCore = false; % C/MR2 is too high to produce a realistic interior structure
% Planet.FeCore = true; % C/MR2 is too high to produce a realistic interior structure

% Seismic.mantleEOS = 'CV3hy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Planet.xFeS_meteoritic = 0.0405; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.55; %0.25
% Planet.xFe_core = 0.0279 ; % this is the total Fe  in Fe and FeS
% Planet.XH2O = 0.0035; % total fraction of water in CV3; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3644;
% Planet.phi_surface = 0;

% Seismic.mantleEOS = 'CIhy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
% Planet.xFeS_meteoritic = 0.0908; %CM2 mean from Jarosewich 1990
% Planet.xFeS = 0.2; %0.25
% Planet.xFe_core = 0.0583 ; % this is the total Fe  in Fe and FeS
% Planet.XH2O = 0.169; % total fraction of water in CI; use this to compute the excess or deficit indicated by the mineralogical model
% Planet.rho_sil_withcore_kgm3 = 3440;


% Seismic.mantleEOS = 'CM2hy1wt_678_1.tab';% (2900 for Q= 100 GW, 3240 for Q= 220 GW)
Planet.xFeS_meteoritic = 0.0676; %CM2 mean from Jarosewich 1990
Planet.xFeS = 1; %0.25
Planet.xFe_core = 0.0463 ; % this is the total Fe  in Fe and FeS
Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model
Planet.rho_sil_withcore_kgm3 = 3000;

% Seismic.mantleEOS = 'CM2hy1wt_678_1.tab';

Params.LineStyle='--';
Params.wrefLine = '--';
% Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [250 255.7]; % 10 Wt% temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);


Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = 255.7; % 10 Wt% temperatures at the bottom of the Ice Ih
outPlanet = PlanetProfile(Planet,Seismic,Params);

asymTables = printInteriorAsym(outPlanet);
fpath_asym = [Planet.name '/interior_model_asym_' Planet.name];
fullPath = [fpath_asym 'Zb' char(outPlanet.ice_thk) 'km.txt'];
writetable(asymTables{1},fullPath);
disp(['Interior conductivity model saved to ' fullPath])

outWaveforms = LayeredInductionResponseJupiter(outPlanet,FTdata,Params);

Params.LineStyle=':';
Planet.Ocean.w_ocean_pct=13; Planet.Tb_K = 257.4; % 10 Wt% temperatures at the bottom of the Ice Ih
outPlanet = PlanetProfile(Planet,Seismic,Params);
outWaveforms = LayeredInductionResponseJupiter(outPlanet,FTdata,Params);

% Params.LineStyle='-';
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [253.1 257.4]; % pure water, temperatures at the bottom of the Ice Ih
% outPlanet = PlanetProfile(Planet,Seismic,Params);

% 
% 
% Planet.FeCore = true;
% Planet.Cmeasured = 0.32; % rounded from the 0.9*.3549; as suggested by Gao and Stevenson 2013
% Planet.Cuncertainty = 0.0042;% Anderson et al. 2001 and Schubert et al. 2004 
% % Seismic.mantleEOS = 'epyro_678_1.tab'; % (3500) this uses the procedure implemented by F. Cammarano
% % Seismic.mantleEOS = 'echon_678_1.tab'; % (3500) this uses the procedure implemented by F. Cammarano
% % Seismic.mantleEOS = 'chon_678_1.tab'; % (3500) this uses the procedure implemented by F. Cammarano
% % Seismic.mantleEOS = 'pyro_678_1.tab'; % (3500) this uses the procedure implemented by F. Cammarano
% % Seismic.mantleEOS = 'chonhp_sat_678.tab'; % (3500) this uses the procedure implemented by F. Cammarano
% % Seismic.mantleEOS = 'echon_hp_sat_PX678_14GPa.tab'; % this uses the procedure implemented by F. Cammarano; this includes Ks and Gs. I had to rerun perlex (6.6.3). not sure why
% % Seismic.mantleEOS = 'pyrohy_678_1.tab'; % (3500) this uses the procedure implemented by F. Cammarano
% % 
% Params.LineStyle='--';
% Params.wrefLine = '--';
% Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = [250 255.7]; % 10 Wt% temperatures at the bottom of the Ice Ih
% Planet.rho_sil_withcore_kgm3 = 3525;
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
% Params.LineStyle='-';
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [253.1 257.4]; % pure water, temperatures at the bottom of the Ice Ih
% Planet.rho_sil_withcore_kgm3 = 3525;
% outPlanet = PlanetProfile(Planet,Seismic,Params);
% % 