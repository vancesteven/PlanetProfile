function PPAriel
%PPAriel
Planet.name='Ariel';
% Planet.rho_kgm3 = 1592; %±150 Jacobson, R. A.; Campbell, J. K.; Taylor, A. H.; Synnott, S. P. (June 1992). "The masses of Uranus and its major satellites from Voyager tracking data and earth-based Uranian satellite data". The Astronomical Journal. 103 (6): 2068?2078. Bibcode:1992AJ....103.2068J. doi:10.1086/116211.
Planet.R_m = 578.9e3; %±0.6 Jacobson, R. A.; Campbell, J. K.; Taylor, A. H.; Synnott, S. P. (June 1992). "The masses of Uranus and its major satellites from Voyager tracking data and earth-based Uranian satellite data". The Astronomical Journal. 103 (6): 2068?2078. Bibcode:1992AJ....103.2068J. doi:10.1086/116211.
Planet.M_kg =1.353e21; %±0.120e21 
Planet.Tsurf_K = 60; 
Planet.Psurf_MPa = 0; 
Planet.Cmeasured = 0.306; %  Hussmann 2006 (not actually measured)
Planet.Cuncertainty = 0.08;% 
Planet.FeCore=false;
    Planet.xFeS = 0.25; %0.25
    Planet.rhoFe = 8000; %8000
    Planet.rhoFeS = 5150; %5150
% Planet.rho_sil_withcore_kgm3 = 2400; % Iess et al. 2014
Planet.rho_sil_withcore_kgm3 = 2700; 

%% &&& Orbital and plotting parameters for use in LayeredInductionResponse
Planet.peaks_Hz = [3.45599e-5 2.30405e-5 1.1519e-5];
Planet.peaks_hr = 1./Planet.peaks_Hz./3600;
Planet.f_orb = 2*pi/2.520/24/3600; % radians per second
Params.wlims = [log(0.001) log(1000)];
% Get Fourier spectrum data
try 
    load(['FTdata' Planet.name],'FTdata');
catch
    dlNow = FTdataMissingWarning();
    if dlNow; load(['FTdata' Planet.name],'FTdata'); end
end
% Planet.ionos_bounds = 100e3;
% Planet.ionosPedersen_sig = 30/100e3;
Planet.ionos_only = [];
Planet.PLOT_SIGS = true;
Planet.ADD_TRANSITION_BOUNDS = false;

% WARNING: The following line was copied from PPCallisto.m because it is
% required by PlanetProfile.m in the current version and does not appear in
% this file. An issue has been opened on GitHub. Delete this comment when
% the value has been corrected.
Planet.XH2O = 0.104; % total fraction of water in CM2; use this to compute the excess or deficit indicated by the mineralogical model


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
Seismic.mantleEOS = 'echon_hp_sat_PX678_14GPa.tab'; % this uses the procedure implemented by F. Cammarano


% Seismic.mantleEOS = 'pyrohp_sat_678_1.tab'; %  (3000) this uses the procedure implemented by F. Cammarano; this includes Ks and Gs. I had to rerun perlex (6.6.3). not sure why
% Seismic.mantleEOSname = 'pyrohpsat';

% Seismic.mantleEOS = 'CV_hhph_DEW17_nofluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CV_hhphD17nofluid';

% Seismic.mantleEOS = 'CV_hhph_DEW17_fluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CV_hhphD17fluid';

% Seismic.mantleEOS = 'CM_hhph_DEW17_nofluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CM_hhphD17nofluid';

% Seismic.mantleEOS = 'CM_hhph_DEW17_fluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CM_hhphD17fluid';

Seismic.mantleEOS = 'CI_hhph_DEW17_nofluid_nomelt_685.tab';
Seismic.mantleEOSname = 'CI_hhphD17nofluid';

% Seismic.mantleEOS = 'CI_hhph_DEW17_fluid_nomelt_685.tab';
% Seismic.mantleEOSname = 'CI_hhphD17fluid';


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
Params.savefigformat = 'epsc';
Params.foursubplots =1;
Params.HOLD = 0; % overlay previous run
Params.Legend = false;
Params.LegendPosition = 'North';
Params.ylim = [910 1170];
Params.Pseafloor_MPa = 30;
Params.nsteps_iceI = 20;
Params.nsteps_ocean = 100; 
Params.nsteps_ref_rho = 30;
Params.nsteps_mantle = 1500;
Params.nsteps_core = 100;
Params.Temps = [245 250 252.5 255 260 265 270];
Params.NOPLOTS = 0; %allows user to limit recreating plots & figures after each run


%% Run the Calculation!
% Planet.Ocean.w_ocean_pct=10;Planet.Tb_K = [272.8 272.9 273 273.1]; % pure water, temperatures at the bottom of the Ice Ih
% PlanetProfile(Planet,Seismic,Params)
% 
% 
% Params.HOLD = true;
% Params.INCLUDE_ELECTRICAL_CONDUCTIVITY = false;
% 
% Params.CALC_NEW=1;
% Params.CALC_NEW_SOUNDSPEEDS=1;
% 
% Planet.Ocean.w_ocean_pct=0; Planet.Tb_K = [273.1 273.15]; % pure water, temperatures at the bottom of the Ice Ih
% PlanetProfile(Planet,Seismic,Params)
% 
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY = 1;
% 
Params.CALC_NEW =1; % Set CALC_NEW options to 0 to re-use profile data when possible. It is recommended to keep CALC_NEW=1 except when intermediate parameters such as layer thicknesses will not change between runs.
Params.CALC_NEW_REFPROFILES=1;
Params.CALC_NEW_SOUNDSPEEDS=1;
% 
Planet.Ocean.comp='MgSO4';
load L_Ice_MgSO4.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');
Params.LineStyle='--';
Params.colororder = 'cm';
Params.wrefLine = '--';
Params.wref=[0 5 10 15];
% 
Planet.Ocean.w_ocean_pct=10; Planet.Tb_K = 269.2;
outPlanet = PlanetProfile(Planet,Seismic,Params);
outWaveforms = LayeredInductionResponse(outPlanet,FTdata,Params);

Planet.Ocean.w_ocean_pct=1; Planet.Tb_K = 273.1;
outPlanet = PlanetProfile(Planet,Seismic,Params);
outWaveforms = LayeredInductionResponse(outPlanet,FTdata,Params);


% running pure water for the MgSO4 case illustrates >1oC error in the Margules parameterization
Params.LineStyle='-';
Params.colororder = 'cm';
Planet.Ocean.w_ocean_pct=0;  Planet.Tb_K =  273.1; % pure water, 
outPlanet = PlanetProfile(Planet,Seismic,Params);
% 
Planet.Ocean.comp='Seawater';
Params.LineStyle='-.';
Params.wref=[0 34 68];
Params.wrefLine = '-.';
Params.colororder = 'cm';

Planet.Ocean.w_ocean_pct=gsw_SSO; Planet.Tb_K = 270.82;
outPlanet = PlanetProfile(Planet,Seismic,Params);
outWaveforms = LayeredInductionResponse(outPlanet,FTdata,Params);

Planet.ALLOW_NEGALPHA = 0;
Planet.Ocean.w_ocean_pct=0.1*gsw_SSO; Planet.Tb_K = 272.57;
outPlanet = PlanetProfile(Planet,Seismic,Params);
outWaveforms = LayeredInductionResponse(outPlanet,FTdata,Params);


Params.HOLD = 1; % overlay previous runs
Params.LineStyle='-';
Params.colororder = 'cm';
Planet.Ocean.w_ocean_pct=0;  Planet.Tb_K =  272.74; % pure water, 
outPlanet = PlanetProfile(Planet,Seismic,Params);

Params.wrefLine =  ':';
Params.wref=[3 5 10];
Planet.Ocean.comp='NH3';
load L_IceNH3_DATA.mat
Planet.Ocean.fnTfreeze_K = griddedInterpolant(PPg',wwg',TT');

Params.CALC_NEW_REFPROFILES=1;
Params.CALC_NEW_SOUNDSPEEDS=1;
Params.HOLD = 1;
Params.LineStyle =  ':';
Planet.Ocean.w_ocean_pct=3; Planet.Tb_K = 269.535; % 0 Wt% temperatures at the bottom of the Ice Ih
%PlanetProfile(Planet,Seismic,Params)