function LayeredInductionResponseJupiter

close all
warning('off','all')
disp('All warnings are turned off. Turn them on again to check for NaN values in the input data.')

 % white background for figures
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultAxesFontSize',16)

CALC_NEW = 0;
DO_EUROPA = 1;
DO_GANYMEDE = 1;
DO_CALLISTO = 1;

PLOT_PERIODS  = 1;
PLOT_WAVEFORMS = 1;
PRINT_TABLES = 1;

%% increments of distance
km = 1e3;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',10*km,'InitialStep',1e-2); % defaults are 1e-3 and 1e-6; maxstep of 10km leads to fast calculations that are satisfactory everwhere to about 0.001%

% w = linspace(log(0.001),log(60),75); used in previous cc
w = linspace(log(0.001),log(1000),75); 
w = exp(w); % cycles per orbit

r0 = km;
y0 = 0;
n = 1;

% previous data
% FieldData = load('GalileanFrequencyData.mat');
% %% switch Bx and By to be consistent with IAU convention used in other work
% deleteme = FieldData.BxFFT_Europa;
% FieldData.BxFFT_Europa = FieldData.ByFFT_Europa;
% FieldData.ByFFT_Europa = deleteme;
% deleteme = FieldData.BxFFT_Ganymede;
% FieldData.BxFFT_Ganymede = FieldData.ByFFT_Ganymede;
% FieldData.ByFFT_Ganymede = deleteme;
% deleteme = FieldData.BxFFT_Callisto;
% FieldData.BxFFT_Callisto = FieldData.ByFFT_Callisto;
% FieldData.ByFFT_Callisto = deleteme;
% clear deleteme;

%% get the FFT data processed by Corey Cochrane. 
% Notes:
% all data are in IAU reference frame of reference

% 10 year duration ...
%    data points: 524288
%    Sep 1, 2018 12:00:00 AM UTC
%    Sep 1, 2028 12:00:00 AM UTC
%    Time Resolution: 0.0016611 Hz, 601.9958 seconds, 10.0333 minutes, 0.16722 hours 
FieldData = ReadGalileanFFTData;
%% set up the peaks to plot and analyze
Europa.peaks_Hz =   [4.946e-5 2.473e-5 3.259e-6];
Ganymede.peaks_Hz = [5.274e-5 2.637e-5 1.617e-6];
Callisto.peaks_Hz = [5.4584e-05 2.7294e-5 6.892e-7];
% Callisto.peaks_Hz = [8.188e-5 2.7294e-5 6.892e-7];

Europa.peaks_hr =   1./Europa.peaks_Hz/3600;
Ganymede.peaks_hr = 1./Ganymede.peaks_Hz/3600;
Callisto.peaks_hr = 1./Callisto.peaks_Hz/3600;
Europa.Name =   'Europa';
Ganymede.Name = 'Ganymede';
Callisto.Name = 'Callisto';

%% create frequency plots similar to those in Seufert et al. 2011
if PLOT_PERIODS
    disp('\begin{table}[]')
    disp('\centering')
    disp('    \hline')
    disp('\begin{tabular}{c | c c c}')
    disp('&     \multicolumn{3}{c}{\textbf{Period (hr)}}\\')
    disp('     &     $B_{x,y,z}$ (nT)&     $B_{x,y,z}$ (nT)&     $B_{x,y,z}$ (nT)\\')
    disp('    \hline')

    nfig = 4000;
    plotPeriods(nfig,Europa,FieldData);
    nfig = 4003;
    plotPeriods(nfig,Ganymede,FieldData);
    nfig = 4005;
    plotPeriods(nfig,Callisto,FieldData);

    disp('   \end{tabular}')
    disp('    \caption{Caption}')
    disp('     \label{tab:my_label}')
    disp('\end{table}')
end

%%
if DO_EUROPA
Re_km = 1561;
r = linspace(km,2*Re_km*1e3,100);

f_orb = 2*pi/3.55/86400; % radians per second
EData.f_orb = f_orb;

    nfig = 125;
    figure(nfig);clf

    %get the field data
EData.frequency = FieldData.frequency;
EData.peaks_Hz = Europa.peaks_Hz;
    
meanEuropa = (FieldData.BxFFT_Europa.^2+FieldData.ByFFT_Europa.^2+FieldData.BzFFT_Europa.^2).^(1/2);
EData.Name = 'Europa';
EData.Bx = FieldData.BxFFT_Europa;
EData.By = FieldData.ByFFT_Europa;
EData.Bz = FieldData.BzFFT_Europa;
EData.meanField_nT = meanEuropa;

if CALC_NEW
    EuropaWaveforms = [];
else
    try
        load(fullfile('Europa','EuropaWaveforms'))
    catch loadWaveformsError
        disp('ERROR: EuropaWaveforms.mat was not found.')
        disp('It probably has not been generated. Set CALC_NEW to 1 to correct this.')
        rethrow(loadWaveformsError)
    end
end

color_warmSW = 	[176,0,255]/255;

%% MgSO4 1wt% 5km ice
if ~isfield(EuropaWaveforms,'M1_5km')
[boundaries,sig,E_M1_5km] = readData('Europa','EuropaProfile_MgSO41WtPct_Ts110Zb4797mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_M1_5km.ionos_bounds = 100e3;
E_M1_5km.ionosPedersen_sig = 30/100e3;
E_M1_5km.ionos_only = [];

E_M1_5km.R = Re_km;
E_M1_5km.PLOT_SIGS = true;
E_M1_5km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 1 wt% - 5 km')
M1_5km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_M1_5km);
M1_5km.Name = '\textbf{MgSO$_{4}$ 1 wt\%}';
EuropaWaveforms.M1_5km = M1_5km;
end
%% MgSO4 1wt% 30 km thick ice
if ~isfield(EuropaWaveforms,'M1_30km')
[boundaries,sig,E_M1_30km] = readData('Europa','EuropaProfile_MgSO41WtPct_Ts110Zb30301mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_M1_30km.ionos_bounds = 100e3;
E_M1_30km.ionosPedersen_sig = 30/100e3;

E_M1_30km.R = Re_km;
E_M1_30km.PLOT_SIGS = false;
E_M1_30km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 1 wt% - 30 km')
M1_30km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_M1_30km);
M1_30km.Name = '';
EuropaWaveforms.M1_30km = M1_30km;
end
%%
%% MgSO4 10wt% 5km ice
if ~isfield(EuropaWaveforms,'M10_5km')
[boundaries,sig,E_M10_5km] = readData('Europa','EuropaProfile_MgSO410WtPct_Ts110Zb4679mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_M10_5km.ionos_bounds = 100e3;
E_M10_5km.ionosPedersen_sig = 30/100e3;

E_M10_5km.R = Re_km;
E_M10_5km.PLOT_SIGS = false;
E_M10_5km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 10 wt% - 5 km')
M10_5km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_M10_5km);
M10_5km.Name = '\textbf{MgSO$_{4}$ 10 wt\%}';
EuropaWaveforms.M10_5km = M10_5km;
end
%% MgSO4 10wt% 30 km ice
if ~isfield(EuropaWaveforms,'M10_30km')
[boundaries,sig,E_M10_30km] = readData('Europa','EuropaProfile_MgSO410WtPct_Ts110Zb29667mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_M10_30km.ionos_bounds = 100e3;
E_M10_30km.ionosPedersen_sig = 30/100e3;

E_M10_30km.R = Re_km;
E_M10_30km.PLOT_SIGS = false;
E_M10_30km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 10 wt% - 30 km')
M10_30km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_M10_30km);
M10_30km.Name = '';
EuropaWaveforms.M10_30km = M10_30km;
end
%% Seawater 3.5ppt 5km ice
if ~isfield(EuropaWaveforms,'Sw_p1_5km')
[boundaries,sig,E_Sw_p1_5km] = readData('Europa','EuropaProfile_Seawater4WtPct_Ts110Zb5085mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_Sw_p1_5km.ionos_bounds = 100e3;
E_Sw_p1_5km.ionosPedersen_sig = 30/100e3;

E_Sw_p1_5km.R = Re_km;
E_Sw_p1_5km.PLOT_SIGS = false;
E_Sw_p1_5km.ADD_TRANSITION_BOUNDS = false; 
disp('==Seawater 0.35165 wt% - 5 km')
Sw_p1_5km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_Sw_p1_5km);
Sw_p1_5km.Name = '\textbf{Seawater 0.35165 wt\%}';
EuropaWaveforms.Sw_p1_5km = Sw_p1_5km;
end
%% Seawater 3.5ppt 30km ice
if ~isfield(EuropaWaveforms,'Sw_p1_30km')
[boundaries,sig,E_Sw_p1_30km] = readData('Europa','EuropaProfile_Seawater4WtPct_Ts110Zb30028mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_Sw_p1_30km.ionos_bounds = 100e3;
E_Sw_p1_30km.ionosPedersen_sig = 30/100e3;

E_Sw_p1_30km.R = Re_km;
E_Sw_p1_30km.PLOT_SIGS = false;
E_Sw_p1_30km.ADD_TRANSITION_BOUNDS = false; 
disp('==Seawater 0.35165 wt% - 30 km')
Sw_p1_30km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_Sw_p1_30km);
Sw_p1_30km.Name = '';
EuropaWaveforms.Sw_p1_30km = Sw_p1_30km;
end
%% Seawater 35ppt 5km ice
if ~isfield(EuropaWaveforms,'Sw_1_5km')
[boundaries,sig,E_Sw_1_5km] = readData('Europa','EuropaProfile_Seawater35WtPct_Ts110Zb4741mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_Sw_1_5km.ionos_bounds = 100e3;
E_Sw_1_5km.ionosPedersen_sig = 30/100e3;

E_Sw_1_5km.R = Re_km;
E_Sw_1_5km.PLOT_SIGS = false;
E_Sw_1_5km.ADD_TRANSITION_BOUNDS = false; 
disp('==Seawater 3.5165 wt% - 5 km')
Sw_1_5km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_Sw_1_5km);
Sw_1_5km.Name = '\textbf{Seawater 3.5165 wt\%}';
EuropaWaveforms.Sw_1_5km = Sw_1_5km;
end
%% Seawater 35ppt 31 km ice
if ~isfield(EuropaWaveforms,'Sw_1_30km')
[boundaries,sig,E_Sw_1_30km] = readData('Europa','EuropaProfile_Seawater35WtPct_Ts110Zb30488mQm6mWm2_CMR2p3460_Europa_CM-CI-67P_silicate_95Fe5S_core_686_fluid_properties_not_included.txt');

E_Sw_1_30km.ionos_bounds = 100e3;
E_Sw_1_30km.ionosPedersen_sig = 30/100e3;

E_Sw_1_30km.R = Re_km;
E_Sw_1_30km.PLOT_SIGS = false;
E_Sw_1_30km.ADD_TRANSITION_BOUNDS = false; 
disp('==Seawater 3.5165 wt% - 5 km')
Sw_1_30km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,EData,E_Sw_1_30km);
Sw_1_30km.Name = '';
EuropaWaveforms.Sw_1_30km = Sw_1_30km;
end
%%
save(fullfile('Europa','EuropaWaveforms'),'EuropaWaveforms')
%% Table output
if PRINT_TABLES
disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
printWaveformTables(EData,EuropaWaveforms,'x');
disp('YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY')
printWaveformTables(EData,EuropaWaveforms,'y');
disp('ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ')
printWaveformTables(EData,EuropaWaveforms,'z');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%% Figures
if PLOT_WAVEFORMS
plotWaveformData(nfig,EData,EuropaWaveforms)
end
%%
% figure(nfig+10);clf;%hold on
% set(gcf,'Position',[157   432   717   583])
% subplot(2,1,1);hold on
% loglog(1./3600./EData.frequency,pdiff(E_M1_5km.InPhaseY,E_M1_5km_mean.InPhaseY),'m--');
% hp(2) = loglog(1./3600./EData.frequency,pdiff(E_M1_30km.InPhaseY,E_M1_30km_mean.InPhaseY),'b--');
% hp(3) = loglog(1./3600./EData.frequency,pdiff(E_M10_5km.InPhaseY,E_M10_5km_mean.InPhaseY),'m--','LineWidth',3);
% hp(4) = loglog(1./3600./EData.frequency,pdiff(E_M10_30km.InPhaseY,E_M10_30km_mean.InPhaseY),'b--','LineWidth',3);
% 
% hp(5) = loglog(1./3600./EData.frequency,pdiff(E_p1Sw_5km.InPhaseY,E_p1Sw_5km_mean.InPhaseY),'c-.');
% hp(6) = loglog(1./3600./EData.frequency,pdiff(E_p1Sw_30km.InPhaseY,E_p1Sw_30km_mean.InPhaseY),'c-.');
% hp(7) = loglog(1./3600./EData.frequency,pdiff(E_Sw_5km.InPhaseY,E_Sw_5km_mean.InPhaseY),'c-.','LineWidth',3);
% hp(8) = loglog(1./3600./EData.frequency,pdiff(E_Sw_31km.InPhaseY,E_Sw_31km_mean.InPhaseY),'c-.','LineWidth',3);
% hvl = vline(1./3600./EData.peaks_Hz,'r');
% set(gca,'xscale','log')
% axis tight
% box on
% % xlabel('Frequency (Hz)')
% ylabel('$\Delta(\mathcal{A}_1^e)_\mathrm{Adiabat,Mean}$ ($\%$)','Interpreter','latex')
% 
% % saveas(gcf,'EuropaMeanVsAdiabat.eps','epsc');
% 
% % figure(nfig+20);clf;hold on
% subplot(2,1,2);hold on
% hp(1) = loglog(1./3600./EData.frequency,pdiff(E_M1_5km.Qinterp,E_M1_5km_top.Qinterp),'m--');
% hp(2) = loglog(1./3600./EData.frequency,pdiff(E_M1_30km.Qinterp,E_M1_30km_top.Qinterp),'b--');
% hp(3) = loglog(1./3600./EData.frequency,pdiff(E_M10_30km.Qinterp,E_M10_30km_top.Qinterp),'b--','LineWidth',3);
% hp(4) = loglog(1./3600./EData.frequency,pdiff(E_M10_5km.Qinterp,E_M10_5km_top.Qinterp),'m--','LineWidth',3);
% 
% hp(5) = loglog(1./3600./EData.frequency,pdiff(E_p1Sw_5km.Qinterp,E_p1Sw_5km_top.Qinterp),'c-.');
% hp(6) = loglog(1./3600./EData.frequency,pdiff(E_p1Sw_30km.Qinterp,E_p1Sw_30km_top.Qinterp),'b-.');
% hp(7) = loglog(1./3600./EData.frequency,pdiff(E_Sw_5km.Qinterp,E_Sw_5km_top.Qinterp),'c-.','LineWidth',3);
% hp(8) = loglog(1./3600./EData.frequency,pdiff(E_Sw_31km.Qinterp,E_Sw_31km_top.Qinterp),'b-.','LineWidth',3);
% 
% hvl = vline(1./3600./EData.peaks_Hz,'r');
% set(gca,'xscale','log')
% xlabel('Period (hr)')
% % xlabel('Frequency (Hz)')
% ylabel('$\Delta(\mathcal{A}_1^e)_\mathrm{Adiabat,Top}$ ($\%$)','Interpreter','latex')
% % ylabel('$(B_\mathrm{Adiabat}-B_\mathrm{Top})/B_\mathrm{Top}$ ($\%$)','Interpreter','latex')
% axis tight
% box on
% % saveas(gcf,'EuropaTopVsAdiabat.eps','epsc');
% saveas(gcf,'EuropaTopMeanVsAdiabat.eps','epsc');

end

if DO_GANYMEDE
    %%%%% Ganymede
Rg_km = 2634;
r = linspace(km,2*Rg_km*1e3,200);

f_orb = 2*pi/7.15/86400; % orbit in Hz
GData.f_orb = f_orb;

    nfig = 126;
    figure(nfig);clf

    %get the mean field data
GData.frequency = FieldData.frequency;
GData.peaks_Hz = Ganymede.peaks_Hz;

meanGanymede = (FieldData.BxFFT_Ganymede.^2+FieldData.ByFFT_Ganymede.^2+FieldData.BzFFT_Ganymede.^2).^(1/2);
GData.Name = 'Ganymede';
GData.Bx = FieldData.BxFFT_Ganymede;
GData.By = FieldData.ByFFT_Ganymede;
GData.Bz = FieldData.BzFFT_Ganymede;
GData.meanField_nT = meanGanymede;

figure(126);clf
%Liu 2016 in situ, measured electrical conductivity of ices: ice V

if CALC_NEW
    GanymedeWaveforms = [];
else
    try
        load(fullfile('Ganymede','GanymedeWaveforms'))
    catch loadWaveformsError
        disp('ERROR: GanymedeWaveforms.mat was not found.')
        disp('It probably has not been generated. Set CALC_NEW to 1 to correct this.')
        rethrow(loadWaveformsError)
    end
end

%create the input from planet profile data
%% MgSO4 1wt% 26km ice
if ~isfield(GanymedeWaveforms,'M1_26km')
[boundaries,sig,G_M1_26km] = readData('Ganymede','GanymedeProfile_MgSO41WtPct_Ts110Zb25614mQm1mWm2_CMR2p3115_CM2hy1wt_678_1.txt');

G_M1_26km.ionos_bounds = 100e3;
G_M1_26km.ionosPedersen_sig = 2/100e3;
G_M1_26km.ionos_only = [];

G_M1_26km.R = Rg_km;
G_M1_26km.PLOT_SIGS = true;
G_M1_26km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 1 wt% - 26 km')
M1_26km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,GData,G_M1_26km);
M1_26km.Name = '\textbf{MgSO$_4$ 1 wt\%}';
GanymedeWaveforms.M1_26km = M1_26km;
end

%% MgSO4 1wt& 95km ice
if ~isfield(GanymedeWaveforms,'M1_95km')
[boundaries,sig,G_M1_95km] = readData('Ganymede','GanymedeProfile_MgSO41WtPct_Ts110Zb94295mQm1mWm2_CMR2p3115_CM2hy1wt_678_1.txt');

G_M1_95km.ionos_bounds = 100e3;
G_M1_95km.ionosPedersen_sig = 2/100e3;

G_M1_95km.R = Rg_km;
G_M1_95km.PLOT_SIGS = true;
G_M1_95km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 1 wt% - 26 km')
M1_95km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,GData,G_M1_95km);
M1_95km.Name = '';
GanymedeWaveforms.M1_95km = M1_95km;
end
%% MgSO4 10wt% 26 km ice
if ~isfield(GanymedeWaveforms,'M10_26km')
[boundaries,sig,G_M10_26km] = readData('Ganymede','GanymedeProfile_MgSO410WtPct_Ts110Zb25240mQm1mWm2_CMR2p3115_CM2hy1wt_678_1.txt');

G_M10_26km.ionos_bounds = 100e3;
G_M10_26km.ionosPedersen_sig = 2/100e3;

G_M10_26km.R = Rg_km;
G_M10_26km.PLOT_SIGS = true;
G_M10_26km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 10 wt% - 26 km')
M10_26km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,GData,G_M10_26km);
M10_26km.Name = '\textbf{MgSO$_4$ 10 wt\%}';
GanymedeWaveforms.M10_26km = M10_26km;
end
%% MgSO4 10wt% 95km ice
% thick ice, the intermediate case from Vance et al. 2018, with 528 km of
% ice VI
if ~isfield(GanymedeWaveforms,'M10_95km')
[boundaries,sig,G_M10_95km] = readData('Ganymede','GanymedeProfile_MgSO410WtPct_Ts110Zb94548mQm1mWm2_CMR2p3115_CM2hy1wt_678_1.txt');

G_M10_95km.ionos_bounds = 100e3;
G_M10_95km.ionosPedersen_sig = 2/100e3;

G_M10_95km.sub_ocean.Depth_km = 500;
G_M10_95km.sub_ocean.Thickness_km = 30;
G_M10_95km.sub_ocean.sig = 20;

G_M10_95km.R = Rg_km;
G_M10_95km.PLOT_SIGS = true;
G_M10_95km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 10 wt% - 95 km')
M10_95km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,GData,G_M10_95km);
M10_95km.Name = '';
GanymedeWaveforms.M10_95km = M10_95km;
end
%%
save(fullfile('Ganymede','GanymedeWaveforms'),'GanymedeWaveforms')
%% add items to the APhi figure
%% Table output
if PRINT_TABLES
disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
printWaveformTables(GData,GanymedeWaveforms,'x');
disp('YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY')
printWaveformTables(GData,GanymedeWaveforms,'y');
disp('ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ')
printWaveformTables(GData,GanymedeWaveforms,'z');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%% Figures
if PLOT_WAVEFORMS
plotWaveformData(nfig,GData,GanymedeWaveforms)
end
%%
% figure(nfig+10);clf;
% set(gcf,'Position',[157   432   717   583])
% subplot(2,1,1);hold on
% hp(1) = loglog(1./3600./GData.frequency,pdiff(G_M1_26km.Qinterp,G_M1_26km_mean.Qinterp),'m--');
% hp(2) = loglog(1./3600./GData.frequency,pdiff(G_M1_95km.Qinterp,G_M1_95km_mean.Qinterp),'b--');
% hp(3) = loglog(1./3600./GData.frequency,pdiff(G_M10_26km.Qinterp,G_M10_26km_mean.Qinterp),'m--','LineWidth',3);
% hp(4) = loglog(1./3600./GData.frequency,pdiff(G_M10_95km.Qinterp,G_M10_95km_mean.Qinterp),'b--','LineWidth',3);
% hvl = vline(1./3600./GData.peaks_Hz,'r');
% set(gca,'xscale','log')
% axis tight
% box on
% % xlabel('Period (hr)')
% % xlabel('Frequency (Hz)')
% ylabel('$\Delta(\mathcal{A}_1^e)_\mathrm{Adiabat,Mean}$ ($\%$)','Interpreter','latex')
% % saveas(gcf,'GanymedeMeanVsAdiabat.eps','epsc');
% 
% 
% % figure(nfig+20);clf;hold on
% subplot(2,1,2);hold on
% hp(1) = loglog(1./3600./GData.frequency,pdiff(G_M1_26km.Qinterp,G_M1_26km_top.Qinterp),'m--');
% hp(2) = loglog(1./3600./GData.frequency,pdiff(G_M1_95km.Qinterp,G_M1_95km_top.Qinterp),'b--');
% hp(3) = loglog(1./3600./GData.frequency,pdiff(G_M10_26km.Qinterp,G_M10_26km_top.Qinterp),'b--','LineWidth',3);
% hp(4) = loglog(1./3600./GData.frequency,pdiff(G_M10_95km.Qinterp,G_M10_95km_top.Qinterp),'m--','LineWidth',3);
% hvl = vline(1./3600./GData.peaks_Hz,'r');
% set(gca,'xscale','log')
% xlabel('Period (hr)')
% % xlabel('Frequency (Hz)')
% ylabel('$\Delta(\mathcal{A}_1^e)_\mathrm{Adiabat,Top}$ ($\%$)','Interpreter','latex')
% axis tight 
% box on
% % saveas(gcf,'GanymedeTopVsAdiabat.eps','epsc');
% saveas(gcf,'GanymedeTopMeanVsAdiabat.eps','epsc');


end

%%%%% Callisto
if DO_CALLISTO
w = linspace(log(0.01),log(100000),100); 
w = exp(w); % cycles per orbit

Rc_km = 2410.3;
r = linspace(km,2*Rc_km*km,100);

f_orb = 2*pi/17/86400; % orbit in Hz
CData.f_orb = f_orb;

nfig = 127;
figure(nfig);clf;

CData.frequency = FieldData.frequency;
indFreqEnd = find(FieldData.frequency>max(Callisto.peaks_Hz));
finds = 1:indFreqEnd(1)+100; % the padding index of 100 is so the max-finder near the end of this file can search within the index space set by the variable hw
CData.frequency = FieldData.frequency(finds);
CData.peaks_Hz = Callisto.peaks_Hz;

    %get the mean field data
    meanCallisto = (FieldData.BxFFT_Callisto.^2+FieldData.ByFFT_Callisto.^2+FieldData.BzFFT_Callisto.^2).^(1/2);
CData.meanField_nT = meanCallisto(finds);
CData.Name = 'Callisto';
CData.Bx = FieldData.BxFFT_Callisto(finds);
CData.By = FieldData.ByFFT_Callisto(finds);
CData.Bz = FieldData.BzFFT_Callisto(finds);

if CALC_NEW
    CallistoWaveforms = [];
else
    try
        load(fullfile('Callisto','CallistoWaveforms'))
    catch loadWaveformsError
        disp('ERROR: CallistoWaveforms.mat was not found.')
        disp('It probably has not been generated. Set CALC_NEW to 1 to correct this.')
        rethrow(loadWaveformsError)
    end
end


if ~isfield(CallistoWaveforms,'M1_100km')
[boundaries,sig,C_M1_100km] = readData('Callisto','CallistoProfile_MgSO41WtPct_Ts110Zb99783mQm2mWm2_CMR2p3549_pyrohp_sat_678_1.txt');

C_M1_100km.ionos_bounds = 100e3;
C_M1_100km.ionosPedersen_sig = 800/100e3;
C_M1_100km.ionosCowling_sig = 6850/100e3;
C_M1_100km.ionos_only = [];

C_M1_100km.R = Rc_km;
C_M1_100km.PLOT_SIGS = true;
C_M1_100km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 1 wt% - 100 km')
M1_100km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,CData,C_M1_100km);
M1_100km.Name = '\textbf{MgSO$_{4}$ 1 wt\%}';
CallistoWaveforms.M1_100km = M1_100km;
end
%% MgSO4 1wt% 130km ice
if ~isfield(CallistoWaveforms,'M1_130km')
[boundaries,sig,C_M1_130km] = readData('Callisto','CallistoProfile_MgSO41WtPct_Ts110Zb129674mQm2mWm2_CMR2p3549_pyrohp_sat_678_1.txt');

C_M1_130km.ionos_bounds = 100e3;
C_M1_130km.ionosPedersen_sig = 800/100e3;
C_M1_130km.ionosCowling_sig = 6850/100e3;

C_M1_130km.R = Rc_km;
C_M1_130km.PLOT_SIGS = true;
C_M1_130km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 1 wt% - 130 km')
M1_130km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,CData,C_M1_130km);
M1_130km.Name = '';
CallistoWaveforms.M1_130km = M1_130km;
end
%% MgSO4 10wt% 100 km ice
if ~isfield(CallistoWaveforms,'M10_100km')
[boundaries,sig,C_M10_100km] = readData('Callisto','CallistoProfile_MgSO410WtPct_Ts110Zb99804mQm2mWm2_CMR2p3549_pyrohp_sat_678_1.txt');

C_M10_100km.ionos_bounds = 100e3;
C_M10_100km.ionosPedersen_sig = 800/100e3;
C_M10_100km.ionosCowling_sig = 6850/100e3;

C_M10_100km.R = Rc_km;
C_M10_100km.PLOT_SIGS = true;
C_M10_100km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 10 wt% - 100 km')
M10_100km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,CData,C_M10_100km);
M10_100km.Name = '\textbf{MgSO$_{4}$ 10 wt\%}';
CallistoWaveforms.M10_100km = M10_100km;
end
%% MgSO4 10wt% 130km ice
if ~isfield(CallistoWaveforms,'M10_130km')
[boundaries,sig,C_M10_130km] = readData('Callisto','CallistoProfile_MgSO410WtPct_Ts110Zb130337mQm2mWm2_CMR2p3549_pyrohp_sat_678_1.txt');

C_M10_130km.ionos_bounds = 100e3;
C_M10_130km.ionosPedersen_sig = 800/100e3;
C_M10_130km.ionosCowling_sig = 6850/100e3;

C_M10_130km.R = Rc_km;
C_M10_130km.PLOT_SIGS = true;
C_M10_130km.ADD_TRANSITION_BOUNDS = false; 
disp('==MgSO4 1 wt% - 130 km')
M10_130km = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,CData,C_M1_130km);
M10_130km.Name = '';
CallistoWaveforms.M10_130km = M10_130km;
end

save(fullfile('Callisto','CallistoWaveforms'),'CallistoWaveforms')
%% add items to the APhi figure
%% Table output
if PRINT_TABLES
disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
printWaveformTables(CData,CallistoWaveforms,'x');
disp('YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY')
printWaveformTables(CData,CallistoWaveforms,'y');
disp('ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ')
printWaveformTables(CData,CallistoWaveforms,'z');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%% Figures
if PLOT_WAVEFORMS
    plotWaveformData(nfig,CData,CallistoWaveforms)
end
%%
% figure(nfig+10);clf;
% set(gcf,'Position',[157   432   717   583])
% subplot(2,1,1);hold on
% hp(1) = loglog(1./3600./CData.frequency,pdiff(C_M1_100km.Qinterp,C_M1_100km_mean.Qinterp),'m--');
% hp(2) = loglog(1./3600./CData.frequency,pdiff(C_M1_130km.Qinterp,C_M1_130km_mean.Qinterp),'b--');
% hp(3) = loglog(1./3600./CData.frequency,pdiff(C_M10_100km.Qinterp,C_M10_100km_mean.Qinterp),'m--','LineWidth',3);
% hp(4) = loglog(1./3600./CData.frequency,pdiff(C_M10_130km.Qinterp,C_M10_130km_mean.Qinterp),'b--','LineWidth',3);
% hvl = vline(1./3600./CData.peaks_Hz,'r');
% set(gca,'xscale','log')
% % xlabel('Period (hr)')
% % xlabel('Frequency (Hz)')
% ylabel('$\Delta(\mathcal{A}_1^e)_\mathrm{Adiabat,Mean}$ ($\%$)','Interpreter','latex')
% axis tight
% box on
% % saveas(gcf,'CallistoMeanVsAdiabat.eps','epsc');
% 
% 
% % figure(nfig+20);clf;hold on
% subplot(2,1,2);hold on
% hp(1) = loglog(1./3600./CData.frequency,pdiff(C_M1_100km.Qinterp,C_M1_100km_top.Qinterp),'m--');
% hp(1) = loglog(1./3600./CData.frequency,pdiff(C_M1_130km.Qinterp,C_M1_130km_top.Qinterp),'b--');
% hp(1) = loglog(1./3600./CData.frequency,pdiff(C_M10_100km.Qinterp,C_M10_100km_top.Qinterp),'b--','LineWidth',3);
% hp(2) = loglog(1./3600./CData.frequency,pdiff(C_M10_130km.Qinterp,C_M10_130km_top.Qinterp),'m--','LineWidth',3);
% hvl = vline(1./3600./CData.peaks_Hz,'r');
% set(gca,'xscale','log')
% xlabel('Period (hr)')
% % xlabel('Frequency (Hz)')
% ylabel('$\Delta(\mathcal{A}_1^e)_\mathrm{Adiabat,Top}$ ($\%$)','Interpreter','latex')
% axis tight
% % set(gca,'xlim',1./3600./[8.6426e-09 3.4885e-05],'ylim',1./3600./[-16.5431 21.8272])
% box on
% % saveas(gcf,'CallistoTopVsAdiabat.eps','epsc');
% saveas(gcf,'CallistoTopMeanVsAdiabat.eps','epsc');

end
warning('on','all');
end % LayeredInductionResponse 

function plotPeriods(nfig,Planet,FData)
    Period = 1./FData.frequency./3600;
    Bx = FData.(['BxFFT_' Planet.Name]);
    By = FData.(['ByFFT_' Planet.Name]);
    Bz = FData.(['BzFFT_' Planet.Name]);

    xlims.Europa = [1.5 2e2];
    xlims.Ganymede = [1 2e3];
    xlims.Callisto = [1 2e3];
   figure(nfig);clf; hold on;    box on
    set(gcf,'Position',[ 73         207        1097         337])
    set(gca,'xscale','log','yscale','log','xlim',xlims.(Planet.Name),'ylim',[1e-5 300])
    text(1.7,40,Planet.Name,'FontSize',40)
    xlabel('Period (hr)');
    ylabel('Component Amplitude (nT)')
    plot(Period,Bx,'b',Period,By,'k',Period,Bz,'g');

    [peaksX_hr,peaksX_nT] = localMax(Period,Bx,Planet.peaks_hr);
    [peaksY_hr,peaksY_nT] = localMax(Period,By,Planet.peaks_hr);
    [peaksZ_hr,peaksZ_nT] = localMax(Period,Bz,Planet.peaks_hr);

    hp = plot(peaksY_hr,peaksY_nT,'ko');
    hp.MarkerFaceColor = 'k';
    for it = 1:length(peaksY_hr)
        ht = text(peaksY_hr(it),peaksY_nT(it),{[num2str(peaksY_hr(it),'%0.2f') ' hr'],['$B_y$: ' num2str(peaksY_nT(it),'%0.2f') ' nT']});
        ht.Interpreter = 'latex';
        ht.FontSize = 16;
        ht.HorizontalAlignment = 'right';
    end

    hl = legend({'$B_x$','$B_y$','$B_z$'},'Interpreter','latex');
    
    % print the tabular information
    dstr = ['\textbf{' Planet.Name '} '];
    for in = 1:3 % count from lowest period (highest frequency)
        dstr = [dstr ' & \textbf{' num2str(peaksY_hr(in),'%0.2f') '}' ];
    end
    disp([dstr ' \\'])
    dstr = '             ';
    for in = 1:3
        dstr = [dstr ' & ' num2str([peaksX_nT(in) peaksY_nT(in) peaksZ_nT(in)],'%8.2f') ];
    end
    disp([dstr ' \\'])
    disp('    \hline')
end %plotPeriods
%%
function [boundaries,sig,Planet] = readData(datpath, savefile)
fdat = importdata(fullfile(datpath,savefile));
nb = 0;
boundaries = [];
sig = [];
% making the wector as small as possible, kiptain
% a small wector is important for solving the diff eq within the lifetime of the
% universe
if strfind(fdat.textdata{end},'phase')
    kind = 12;
else
    kind = 10;
end
for ij = length(fdat.data(:,1))-1:-1:1
    if isnan(fdat.data(ij+1,kind))
        warning(['found NaN in conductivity data read from file at entry ' num2str(ij+1) ' of ' num2str(length(fdat.data(:,1)-1)) '.'])
    elseif fdat.data(ij+1,kind)~=fdat.data(ij,kind) %
        nb = nb+1;
        boundaries(nb) = fdat.data(ij+1,3);
        sig(nb) = fdat.data(ij+1,kind);
    end
end
% upper boundary
boundaries(nb+1) = fdat.data(1,3);
sig(nb+1) = fdat.data(1,kind);

sig(sig==0)=1e-16; % zeros make the integration unstable
sig(isnan(sig))=1e-16;
boundaries = boundaries*1000; % in km

% now get the planetary properties to put into a table for publication
% search for liquid where density is greater than 1000 kg/m3
inds = find(fdat.data(:,4)>1e3);
% ld = length(fdat.data(:,4)); % I don't think this is used.

inds2 = find(diff(fdat.data(:,4))>10);
Planet.Tb_K = fdat.data(inds(2),2);
Planet.Tmean_K = mean(fdat.data(inds(1):inds2(2),2));
if strfind(fdat.textdata{end},'phase')
    phase = fdat.data(2:end,11);
    inds_iceIh = find(phase==1);
    inds_ocean = find(phase==0);
    Planet.D_Ih_km = fdat.data(1,3)-fdat.data(inds_iceIh(end),3);
    Planet.D_ocean_km = fdat.data(inds_ocean(1),3)-fdat.data(inds_ocean(end),3);    
else
    Planet.D_Ih_km = fdat.data(1,3)-fdat.data(inds(1),3);
    Planet.D_ocean_km = fdat.data(inds(1),3)-fdat.data(inds2(2),3);
end
Planet.kmean = mean(fdat.data(inds(1):inds2(2),kind));

end % readData
%%
function outPlanet = getSw(r,n,w,sig,boundaries,R_m,r0,y0,opts,FData,Planet)
    opts.hiprec = true;
    opts.par = false;
%     Q = MagComplexA(boundaries,sig,nE*w,n,opts);
    [k,Q] = getMagResponseFunction(r,n,FData.f_orb,w,sig,boundaries,R_m,r0,y0,opts);
    nbound = length(boundaries(:,1));
    outPlanet = Planet;
        for lm = 1:nbound
            Amplitude(lm,:) = 2*abs(squeeze(Q(:,:,lm)));
            if isfield(Planet,'SURF_NORM') && Planet.SURF_NORM
                Amplitude = Amplitude*(Planet.Rtop/Planet.R)^3;
            end
            Phase(lm,:) = -angle(squeeze(Q(:,:,lm)));                
        end
        outPlanet.Amplitude = Amplitude;
        outPlanet.Phase = Phase;
    outPlanet.Q = Q;
    outPlanet.k = k;
    outPlanet.w = w;
    outPlanet.sig = sig;
    outPlanet.boundaries = boundaries;
end % getSw    
%%
function outPlanet = getPlanetMagWaveforms(r,n,w,sig,boundaries,r0,y0,opts,FData,Planet)
% get Sw for the different configurations related to the main input model
if Planet.ADD_TRANSITION_BOUNDS
    % Adding discrete transitions makes a more intuitive looking plot of the sigma values. 
    % Setting the transition values to eps seems to change the output InPhase and Quad values by about 0.1%, while setting the value to r0 
    % doesn't seem to change them at all. I am assuming that increments smaller than step size r0 are introducing errors in the ode solver
    btrans = r0; % or eps
    boundaries = [boundaries(1) boundaries(2)-btrans boundaries(2:end-1) boundaries(end-1)+btrans boundaries(end)];  
    sig = [sig(1) eps sig(2:end-1) eps sig(end)];

    b_mionos = [boundaries boundaries(end)+btrans Planet.R*1e3+Planet.ionos_bounds];
    sPedersen_mionos = [sig Planet.ionosPedersen_sig(1) Planet.ionosPedersen_sig];

    b_ionos = [boundaries boundaries(end)+btrans Planet.R*1e3+Planet.ionos_bounds];
    sPedersen_ionos = [eps*ones(1,length(sig))  Planet.ionosPedersen_sig(1) Planet.ionosPedersen_sig];
else
    b_mionos = [boundaries Planet.R*1e3+Planet.ionos_bounds];
    sPedersen_mionos = [sig Planet.ionosPedersen_sig];

    b_ionos = [boundaries Planet.R*1e3+Planet.ionos_bounds];
    sPedersen_ionos = [1e-16*ones(1,length(sig)) Planet.ionosPedersen_sig];
end

if isfield(Planet,'sub_ocean')
    Dice = Planet.sub_ocean.Depth_km*1e3;
    Dlayer = Planet.sub_ocean.Thickness_km*1e3;
    sig_layer = Planet.sub_ocean.sig;
    b_subo = [boundaries(1)-Dice-Dlayer boundaries(1)-Dice boundaries];
    sPedersen_subo = [1e-16,sig_layer,1e-16,sig(2:end-1),1e-16];
    
    b_subo_ionos = [b_subo Planet.R*1e3+Planet.ionos_bounds];
    sPedersen_subo_ionos = [sPedersen_subo Planet.ionosPedersen_sig];
end

if isfield(Planet,'ionosCowling_sig')
    sCowling_mionos = [sig Planet.ionosCowling_sig];
    sCowling_ionos = [1e-16*ones(1,length(sig)) Planet.ionosCowling_sig];
end
    ocean_sigs = sig(sig>eps);
    [s_mean,s_top] = deal(sig);

    % s_mean = [eps,[1 1]*mean(sig(2:end-1)),eps];
    % b_mean = boundaries([1,round(length(boundaries(2:end-1))/2),end-1,end]);
    b_mean = boundaries;
    s_mean(sig>eps) = mean(ocean_sigs);

    % s_top = [eps,[1 1]*sig(end-1),eps];
    % b_top = boundaries([1,round(length(boundaries(2:end-1))/2),end-1,end]);
    b_top = boundaries;
    s_top(sig>eps) = ocean_sigs(end);


if Planet.PLOT_SIGS
    figure (3000);clf; box on; hold on;
    plot(sig,boundaries/1e3,...
        sPedersen_mionos,b_mionos/1e3,...
        sPedersen_ionos,b_ionos/1e3,...
        s_mean,b_mean/1e3,...
        s_top,b_top/1e3);
    
    if isfield(Planet,'sub_ocean')
        plot(sPedersen_subo,b_subo/1e3,'--',...
            sPedersen_subo_ionos,b_subo_ionos,'--');
    end
    legend({'main field','with ionosphere','ionosphere only','mean ocean','top ocean'},'NumColumns',2,'Location','northwest');
%     legend('boxoff')
    xlabel('\sigma (S/m)')
    ylabel('r (km)')
    set(gca,'xscale','log')
    xlim([1e-6 Inf])
end

Planet.SURF_NORM = false;
disp(' (o)  getting main field');outPlanet.main = getSw(r,n,w,sig,boundaries,Planet.R*1e3,r0,y0,opts,FData,Planet);

Planet.SURF_NORM = true; 
Planet.Rtop = b_mionos(end)/1e3;
disp(' +o+  with Pedersen ionsphere');outPlanet.ionosPedersen = getSw(r,n,w,sPedersen_mionos,b_mionos,b_mionos(end),r0,y0,opts,FData,Planet);
if isfield(Planet,'ionosCowling_sig')
    disp('-+o+- with Cowling ionsphere');outPlanet.ionosCowling = getSw(r,n,w,sCowling_mionos,b_mionos,b_mionos(end),r0,y0,opts,FData,Planet);
end
if isfield(Planet,'ionos_only') % only create the ionosphere without the ocean in the few places it's specified among the inputs
    Planet.Rtop = b_mionos(end)/1e3;
    disp(' + +  Pedersen ionsphere only');    outPlanet.ionos_onlyPedersen = getSw(r,n,w,sPedersen_ionos,b_ionos,b_ionos(end),r0,y0,opts,FData,Planet);
    if isfield(Planet,'ionosCowling_sig')
        disp('-+o+- Cowling ionsphere only');outPlanet.ionos_onlyCowling = getSw(r,n,w,sCowling_ionos,b_mionos,b_mionos(end),r0,y0,opts,FData,Planet);
    end
end


Planet.SURF_NORM = false;
disp('  -   mean ocean');outPlanet.mean = getSw(r,n,w,s_mean,b_mean,Planet.R*1e3,r0,y0,opts,FData,Planet);
disp('  ^   top ocean');outPlanet.top = getSw(r,n,w,s_top,b_top,Planet.R*1e3,r0,y0,opts,FData,Planet);

if isfield(Planet,'sub_ocean')
    disp(' - -  two oceans');    outPlanet.sub_ocean = getSw(r,n,w,sPedersen_subo,b_subo,Planet.R*1e3,r0,y0,opts,FData,Planet);
    Planet.SURF_NORM = true; 
    disp(' -+-  two oceans, with ionosphere');    outPlanet.sub_ocean_ionos = getSw(r,n,w,sPedersen_subo_ionos,b_subo_ionos,b_subo_ionos(end),r0,y0,opts,FData,Planet);    
end
end %getPlanetMagWaveforms

%% Output the data as LaTeX tables. The \end{} parts of the tables need to be supplied
function printWaveformTables(FData,Waveforms,dir)
makeTableHeader(FData,dir)
fnames = fieldnames(Waveforms);
for in = 1:length(fnames)
    plot_opts = getPlotOpts(FData.Name,fnames{in});
    if isfield(plot_opts,'Tb_colorTex')
        table_opts.Tb_colorTex = plot_opts.Tb_colorTex;
    end
    table_opts.dir = dir;
    table_opts.main = Waveforms.(fnames{in}).main; % used for printing PERCENT values
    if isfield(Waveforms.(fnames{in}),'ionos_onlyPedersen')
        table_opts.IONOSPHERE_ONLY = true;
        table_opts.INCLUDE_TDs = false;
        printTableStr('Pedersen',FData,Waveforms.(fnames{in}).ionos_onlyPedersen,table_opts);    
        if isfield(Waveforms.(fnames{in}),'ionos_onlyCowling')
            table_opts.IONOSPHERE_ONLY = false;
            printTableStr('Cowling',FData,Waveforms.(fnames{in}).ionos_onlyCowling,table_opts);                
        end
        disp(['\hline']) % end the section providing ionosphere strengths
    end
    table_opts.IONOSPHERE_ONLY = false;

    % adiabat--this is the main data source. others can be listed in
    % percent of this value
    table_opts.INCLUDE_TDs = true;
    printTableStr(Waveforms.(fnames{in}).Name,FData,Waveforms.(fnames{in}).main,table_opts);
    table_opts.INCLUDE_TDs = false;

    if strcmp(FData.Name,'Callisto')
      table_opts.PERCENT = 0;
    else
      table_opts.PERCENT = 1;
    end
        % list the results for the main calculation, adding the ionospheric
    % effects
    printTableStr('Pedersen',FData,Waveforms.(fnames{in}).ionosPedersen,table_opts);
    if isfield(Waveforms.(fnames{in}),'ionosCowling')
        printTableStr('Cowling',FData,Waveforms.(fnames{in}).ionosCowling,table_opts);                
    end
    
    table_opts.PERCENT = 1;
    %mean ocean
    table_opts.main = Waveforms.(fnames{in}).main;
    printTableStr(['$\overline{\sigma}=$ ' num2str(Waveforms.(fnames{in}).mean.sig(end-3),'%0.4f') ' S/m'],FData,Waveforms.(fnames{in}).mean,table_opts);
    %top ocean
    printTableStr(['$\sigma_\mathrm{top}=$ ' num2str(Waveforms.(fnames{in}).top.sig(end-3),'%0.4f') ' S/m'],FData,Waveforms.(fnames{in}).top,table_opts);
    table_opts.PERCENT = 0;

    table_opts.PERCENT = 1;
    % any modeled second ocean under high-pressure ice
    if isfield(Waveforms.(fnames{in}),'sub_ocean')
        d_layer = Waveforms.(fnames{in}).sub_ocean.sub_ocean.Thickness_km;
        sig_layer = Waveforms.(fnames{in}).sub_ocean.sub_ocean.sig;
         printTableStr(['\textbf{bottom layer: ' num2str(d_layer) ' km ' num2str(sig_layer) ' S/m}'],FData,Waveforms.(fnames{in}).sub_ocean,table_opts); 
         printTableStr('Pedersen',FData,Waveforms.(fnames{in}).sub_ocean_ionos,table_opts); 
    end
    table_opts.PERCENT = 0;
    disp('\hline')
end
end %printWaveformTables
%%
function plotWaveformData(nfig,FData,Waveforms)
fnames = fieldnames(Waveforms);
for in = 1:length(fnames)
    plot_opts = getPlotOpts(FData.Name,fnames{in});
    plotAPhi(nfig,plot_opts,FData,Waveforms.(fnames{in}));
end

% Add labels and save the Aphi figure
figure(nfig*20);
labelWaveformFigs('x',FData.Name)
figure(nfig*20+1);
labelWaveformFigs('y',FData.Name)
figure(nfig*20+2);
labelWaveformFigs('z',FData.Name)
end %plotWaveformData
%% Set the line and marker styles. Conditional formatting strongly depends on the choice and labeling of model inputs. This needs to be generalized for future work.
function plot_opts = getPlotOpts(pname,fname)
switch pname
    case 'Europa'
        switch fname(1) 
        case 'M'
            plot_opts.line='--';
            if strfind(fname,'30')
                [plot_opts.LC,plot_opts.MEC] = deal('b');
                plot_opts.point = '^';
                plot_opts.Tb_colorTex = '\color{blue}';
            else
                [plot_opts.LC,plot_opts.MEC] = deal('m');
                plot_opts.point = 'v';
                plot_opts.Tb_colorTex = '\color{magenta}';
            end
            if strfind(fname,'10')
                plot_opts.MFC = plot_opts.MEC;
                plot_opts.MEC = 'k';
                plot_opts.LW = 3;
            else
                plot_opts.MFC = 'none';
                plot_opts.LW = 1;
            end
        case 'S'
            color_warmSW = 	[176,0,255]/255;
            plot_opts.line='-';
            if strfind(fname,'30')
                [plot_opts.LC,plot_opts.MEC] = deal('c');
                plot_opts.point = '^';
                plot_opts.Tb_colorTex = '\color{cyan}';
            else
               [plot_opts.LC,plot_opts.MEC] = deal(color_warmSW);
               plot_opts.point = 'v';
               plot_opts.Tb_colorTex = '\color[HTML]{b000ff}';
            end
            if strfind(fname,'_1')
                plot_opts.MFC = plot_opts.MEC;
                plot_opts.MEC = 'k';
                plot_opts.LW = 3;
            else
                plot_opts.MFC = 'none';
                plot_opts.LW = 1;
            end
        end
    case 'Ganymede'
        plot_opts.line='--';
        if strfind(fname,'95')
            [plot_opts.LC,plot_opts.MEC] = deal('b');
            plot_opts.point = '^';
            plot_opts.Tb_colorTex = '\color{blue}';
        else
            [plot_opts.LC,plot_opts.MEC] = deal('m');
            plot_opts.point = 'v';
            plot_opts.Tb_colorTex = '\color{magenta}';
        end
        if strfind(fname,'10')
            plot_opts.MFC = plot_opts.MEC;
            plot_opts.MEC = 'k';
            plot_opts.LW = 3;
        else
            plot_opts.MFC = 'none';
            plot_opts.LW = 1;
        end
    case 'Callisto'
        plot_opts.line='--';
        if strfind(fname,'130')
            [plot_opts.LC,plot_opts.MEC] = deal('b');
            plot_opts.point = '^';
            plot_opts.Tb_colorTex = '\color{blue}';
        else
            [plot_opts.LC,plot_opts.MEC] = deal('m');
            plot_opts.point = 'v';
            plot_opts.Tb_colorTex = '\color{magenta}';
        end
        if strfind(fname,'10_')
            plot_opts.MFC = plot_opts.MEC;
            plot_opts.MEC = 'k';
            plot_opts.LW = 3;
        else
            plot_opts.MFC = 'none';
            plot_opts.LW = 1;
        end
end
end %getPlotOpts
%%
function plotDiagnostics(nfig,r,w,PName)
xscale = 'linear';
figure(50+nfig);
subplot(2,2,1)
surf(r/1e3,w,real(k))
set(gca,'xscale',xscale,'zscale','log')
axis tight; 
ylabel(['frequency (cycles per ' PName ' orbit)']);
xlabel('r (km)')
zlabel('Re(k)')

subplot(2,2,2)
surf(r/1e3,w,imag(k));
set(gca,'xscale',xscale,'zscale','log')
axis tight
ylabel(['frequency (cycles per ' PName ' orbit)']);
xlabel('r (km)')
zlabel('Im(k)')

subplot(2,2,3)
surf(r/1e3,w,real(k))
set(gca,'xlim',[1400 1561],'zscale','log')
ylabel(['frequency (cycles per ' PName ' orbit)']);
xlabel('r (km)')
zlabel('Re(k)')
% axis tight
subplot(2,2,4)
surf(r/1e3,w,imag(k));
set(gca,'xlim',[1400 1561],'zscale','log')
ylabel(['frequency (cycles per ' PName ' orbit)']);
xlabel('r (km)')
zlabel('Im(k)')    
end %plotDiagnostics

%% main plots ++++++++++
function plotAPhi(nfig,plot_opts,FData,Planet)
%     global FREQUENCY_AXIS
    f_orb = FData.f_orb;
    f_small = f_orb/2/pi*Planet.main.w; % frequency spectrum
    p_small = 1./f_small./3600;
    Period = 1./FData.frequency./3600;
    PeriodsOfInterest = 1./FData.peaks_Hz./3600;
    
    LineWidth = plot_opts.LW;

    itype = 'spline';
    Amp = interp1(p_small,Planet.main.Amplitude,Period,itype);
    Phase = interp1(p_small,Planet.main.Phase,Period,itype);
    A1e_adiabat = interp1(p_small,Planet.main.Q,Period,itype);
    A1e_mean = interp1(p_small,Planet.mean.Q,Period,itype);
    A1e_top = interp1(p_small,Planet.top.Q,Period,itype);

     ReA1e = Amp.*cos(Phase);
     ImA1e = Amp.*sin(Phase);
%     ReA1e = real(A1e_adiabat);
%     ImA1e = imag(A1e_adiabat);

    InPhaseX = FData.Bx.*ReA1e;
    InPhaseY = FData.By.*ReA1e;
    InPhaseZ = FData.Bz.*ReA1e;

    QuadX = FData.Bx.*ImA1e;
    QuadY = FData.By.*ImA1e;
    QuadZ = FData.Bz.*ImA1e;

    figure(nfig);
    subplot(2,1,1);hold on; hl = [];
            hthis = semilogx(Period,Amp,'LineWidth',LineWidth,'LineStyle',plot_opts.line,'Color',plot_opts.LC); 
            xlabel(['Period (hr)']);
        ylabel('Normalized Amplitude')
        axis tight
    
    ylim([0 1])
    yticks(0:0.2:1)
            hv = vline(PeriodsOfInterest,'r'); 
%             hv = vline(nE/2/pi*wOfInterest,'r'); 
        [hv.Color] = deal('r');
        [hv.LineWidth] = deal(1);
        set(gca,'xscale','log')
    box on

    subplot(2,1,2); hold on
            hpf = semilogx(Period,Phase*180/pi,'LineWidth',LineWidth,'LineStyle',plot_opts.line,'Color',plot_opts.LC);
            xlabel(['Period (hr)']);
        ylabel('Phase Delay (degrees)')
        axis tight
    ylim([0 90])
            hv = vline(PeriodsOfInterest,'r'); 
        [hv.LineWidth] = deal(1);
        set(gca,'xscale','log')
    box on    

    hold on
    switch FData.Name % set the zoom level for each plot
        case 'Europa'
            plot_opts.Bxxlim = [0 12];
            plot_opts.Bxylim = [0 1.7];
            plot_opts.Byxlim = [0 16];
            plot_opts.Byylim = [0 5];
            plot_opts.Bzxlim = [3 10];
            plot_opts.Bzylim = [0 0.7];
        case 'Ganymede'
            plot_opts.Bxxlim = [0 1.8];
            plot_opts.Bxylim = [0 0.12];
            plot_opts.Byxlim = [0 2.5];
            plot_opts.Byylim = [0 0.6];
            plot_opts.Bzxlim = [0 1.8];
            plot_opts.Bzylim = [0 0.08];
        case 'Callisto'
            plot_opts.Bxxlim = [0 0.008];
            plot_opts.Bxylim = [0 0.01];
            plot_opts.Byxlim = [0 0.2];
            plot_opts.Byylim = [0 0.3];
            plot_opts.Bzxlim = [0 0.002];
            plot_opts.Bzylim = [0 0.012];
    end
        
    if ~isfield(plot_opts,'MarkerSize')
        plot_opts.MarkerSize = 10;
    end
        [~,peaksInPhaseX_nT] = localMax(Period,InPhaseX,PeriodsOfInterest);
        [~,peaksInPhaseY_nT] = localMax(Period,InPhaseY,PeriodsOfInterest);
        [~,peaksInPhaseZ_nT] = localMax(Period,InPhaseZ,PeriodsOfInterest);

        [~,peaksQuadX_nT] = localMax(Period,QuadX,PeriodsOfInterest);
        [~,peaksQuadY_nT] = localMax(Period,QuadY,PeriodsOfInterest);
        [~,peaksQuadZ_nT] = localMax(Period,QuadZ,PeriodsOfInterest);
        plot_opts.MS = [4 7 12];
        plot_opts.symbols = '  o';

        plotInductionWaveforms(peaksInPhaseX_nT,peaksQuadX_nT,plot_opts,'x',nfig*20);
        plotInductionWaveforms(peaksInPhaseY_nT,peaksQuadY_nT,plot_opts,'y',nfig*20+1);
        plotInductionWaveforms(peaksInPhaseZ_nT,peaksQuadZ_nT,plot_opts,'z',nfig*20+2);
end %plotAPhi
%%
function plotInductionWaveforms(pInPhase,pQuad,plot_opts,dir,nfig)
figure(nfig);set(gcf,'Position',[469   767   589   251])
subplot(1,2,1);hold on
for ip = 3 % plot circles underneath the orbital period
    hr = plot(pInPhase(ip),pQuad(ip),plot_opts.symbols(ip),'MarkerFaceColor',plot_opts.MFC,'MarkerEdgeColor',plot_opts.MEC,'MarkerSize',plot_opts.MS(ip));
end
for ip = 1:length(pInPhase) % this is intended to break if more than three frequencies are specified
    hr = plot(pInPhase(ip),pQuad(ip),plot_opts.point,'MarkerFaceColor',plot_opts.MFC,'MarkerEdgeColor',plot_opts.MEC,'MarkerSize',plot_opts.MS(ip));
end
xlim(plot_opts.(['B' dir 'xlim']));
ylim(plot_opts.(['B' dir 'ylim']));
xlabel(['Re (nT)'],'Interpreter','latex')
ylabel(['Im (nT)'],'Interpreter','latex')
box on; grid on

subplot(1,2,2);hold on
for ip = 3 % plot circles underneath the orbital period
    hr = plot(pInPhase(ip),pQuad(ip),plot_opts.symbols(ip),'MarkerFaceColor',plot_opts.MFC,'MarkerEdgeColor',plot_opts.MEC,'MarkerSize',plot_opts.MS(ip));
end
for ip = 1:length(pInPhase) % this is intended to break if more than three frequencies are specified
    hr = plot(pInPhase(ip),pQuad(ip),plot_opts.point,'MarkerFaceColor',plot_opts.MFC,'MarkerEdgeColor',plot_opts.MEC,'MarkerSize',plot_opts.MS(ip));
end
set(gca,'YAxisLocation','right')
xlabel(['Re (nT)'],'Interpreter','latex')
ylabel(['Im (nT)'],'Interpreter','latex')
axis tight; box on; grid on
end %plotInductionWaveforms
%% puts labels on the plots of the real and imaginary components of the induction response
function labelWaveformFigs(dir,Name)
hf = gcf;
switch Name
    case 'Europa'
    gfxdir = 'Europa/figures/';
    ht = title(['$|B_{' dir '}|\mathcal{A}_1^e$, $\mathcal{A}_1^e=Ae^{-i\phi}$'],'Interpreter','latex','FontSize',18);
    % [ -6.4312   12.4896         0]);

    annotation(hf,'textbox',...
        [0.143509135200974 0.717573221757318 0.222289890377588 0.190376569037657],...
        'String',{'Europa'},...
        'FontSize',30,...
        'EdgeColor','none');
    % Create textarrows
    annotation(hf,'textarrow',[0.579780755176614 0.595615103532278],...
        [0.594142259414226 0.527196652719665],'String',{'synodic'},'FontSize',15);
    annotation(hf,'textarrow',[0.736906211936663 0.721071863580999],...
        [0.757322175732217 0.694560669456067],'String',{'orbital'},'FontSize',15);
    annotation(hf,'textarrow',[0.87088915956151 0.892813641900122],...
        [0.744769874476987 0.820083682008368],'String',{'synodic','harmonic'},...
        'FontSize',15);
    saveas(hf,fullfile([gfxdir 'AeEuropa' dir '.fig']));
    case 'Ganymede'
        gfxdir = 'Ganymede/figures/';
        annotation(hf,'textbox',...
        [0.143509135200974 0.717573221757318 0.222289890377588 0.190376569037657],...
        'String',{'Ganymede'},...
        'FontSize',30,...
        'EdgeColor','none');
    % Create textarrows
    annotation(hf,'textarrow',[0.579780755176614 0.595615103532278],...
        [0.594142259414226 0.527196652719665],'String',{'synodic'},'FontSize',15);
    annotation(hf,'textarrow',[0.736906211936663 0.721071863580999],...
        [0.757322175732217 0.694560669456067],'String',{'orbital'},'FontSize',15);
    annotation(hf,'textarrow',[0.87088915956151 0.892813641900122],...
        [0.744769874476987 0.820083682008368],'String',{'synodic','harmonic'},...
        'FontSize',15);
    saveas(hf,fullfile([gfxdir 'AeGanymede' dir '.eps']),'epsc');
    case 'Callisto'
        gfxdir = 'Callisto/figures/';
        annotation(hf,'textbox',...
        [0.143509135200974 0.717573221757318 0.222289890377588 0.190376569037657],...
        'String',{'Callisto'},...
        'FontSize',30,...
        'EdgeColor','none');
    % Create textarrows
    annotation(hf,'textarrow',[0.579780755176614 0.595615103532278],...
        [0.594142259414226 0.527196652719665],'String',{'synodic'},'FontSize',15);
    annotation(hf,'textarrow',[0.736906211936663 0.721071863580999],...
        [0.757322175732217 0.694560669456067],'String',{'orbital'},'FontSize',15);
    annotation(hf,'textarrow',[0.87088915956151 0.892813641900122],...
        [0.744769874476987 0.820083682008368],'String',{'synodic','harmonic'},...
        'FontSize',15);
    saveas(hf,fullfile([gfxdir 'AeCallisto' dir '.eps']),'epsc');
end
end %labelWaveformFigs
%% formatting for output
function makeTableHeader(Data,dir)
disp('\begin{table}[h]');
disp('\begin{center}');

d_str = '\begin{tabular}{c c c c|';
% d_str = '\begin{tabular}{|c|c|c|c|c|';
for iF = 1:length(Data.peaks_Hz)
    d_str = [d_str 'c '];
end
    d_str = [d_str '}'];
disp(d_str);

disp('\hline');

printTableFreqs(Data,dir)

n_peaks = length(Data.peaks_Hz);
disp([' $T_{b}$ & $\overline{T}$ &  $D_\mathrm{I}$  & $D_\mathrm{ocean}$        &  \multicolumn{' num2str(n_peaks) '}{c}{$B_' dir '\mathcal{A}_1^e$} \\']);
disp(['  (K)     & (K)       &    (km)     & (km)               & \multicolumn{' num2str(n_peaks) '}{c}{(nT)} \\']);
disp('\hline');
end %makeTableHeader
%%
function printTableFreqs(Data,dir)
disp('\hline');
amps = [];
for iF = 1:length(Data.peaks_Hz)-1 % this will fail if f_inds has only one entry
    amps = [amps ' &'];
end
    amps = [amps ' \\'];

    
[PeakPeriods_hr,PeakFields_nT] = localMax(1./Data.frequency/3600,Data.(['B' dir]),1./Data.peaks_Hz/3600);
    
%print the peak frequencies
d_str = [ '\textbf{' Data.Name '} &   &  &   \multicolumn{1}{r}{Period (hr):}  '];
%  d_str = [Name ' &         &          &                         \multicolumn{2}{c}{$f$ ($\times 10^{-6}$Hz)}  '];
for iF = 1:length(Data.peaks_Hz)
    d_str = [d_str ' & \textbf{' num2str(PeakPeriods_hr(iF),'%0.2f') '}'];
%     d_str = [d_str ' & \textbf{' num2str(1e6*(Data.peaks_Hz(iF)),'%0.2f') '}'];
end
disp([d_str ' \\']);

% now print the associated magnitudes of the peak field
 d_str = [ ' &   &  &   \multicolumn{1}{r}{$B_' dir '$ (nT):}  '];
%  d_str = [Name ' &         &          &                         \multicolumn{2}{c}{$f$ ($\times 10^{-6}$Hz)}  '];
for iF = 1:length(Data.peaks_Hz)
    d_str = [d_str ' & \textbf{' num2str(PeakFields_nT(iF),'%0.2f') '}'];
%     d_str = [d_str ' & \textbf{' num2str(1e6*(Data.peaks_Hz(iF)),'%0.2f') '}'];
end
disp([d_str ' \\']);
disp('\hline');
end % printTableFreqs
%%
function printTableStr(Name,FData,Planet,opts)
    f_orb = FData.f_orb;
    f_small = f_orb/2/pi*Planet.w;
    Period = 1./FData.frequency./3600;
    PeakPeriods = 1./FData.peaks_Hz/3600;
    
    Amp = interp1(f_small,Planet.Amplitude,FData.frequency,'spline');
    Phase = interp1(f_small,Planet.Phase,FData.frequency,'spline');

    if isfield(opts,'PERCENT') && opts.PERCENT ==1
        format = '%0.2f'; % output to 0.01% precision
        Amp_main = interp1(f_small,opts.main.Amplitude,FData.frequency,'spline');
        Phase_main = interp1(f_small,opts.main.Phase,FData.frequency,'spline');
        InPhase = pdiff(Amp.*cos(Phase),Amp_main.*cos(Phase_main));
        Quad = pdiff(Amp.*sin(Phase),Amp_main.*sin(Phase_main));
    else
        format = '%0.3f';
        InPhase = Amp.*cos(Phase).*FData.(['B' (opts.dir)]);
        Quad = Amp.*sin(Phase).*FData.(['B' (opts.dir)]);  
    end
    
if opts.IONOSPHERE_ONLY    
        disp(['\multicolumn{4}{l}{\textbf{Ionosphere Only}} & \textbf{Re Im} & \textbf{Re Im} & \textbf{Re Im}\\'])
        disp('\hline')
end   
if opts.INCLUDE_TDs
    if ~isempty(Name)
        disp(['\multicolumn{4}{l}{' Name '} & \textbf{Re Im} & \textbf{Re Im} & \textbf{Re Im}\\'])
        disp(['\hline'])
    end
    if isfield(opts,'Tb_colorTex')
        col = opts.Tb_colorTex;
    else
        col = '';
    end
    d_str = ['{' col num2str(Planet.Tb_K,'%0.1f') '}'];
    d_str = [d_str ' &' num2str(Planet.Tmean_K,'%0.1f')];
    d_str = [d_str ' &' num2str(Planet.D_Ih_km,'%0.0f')];
    d_str = [d_str ' &' num2str(Planet.D_ocean_km,'%0.0f')];
elseif isfield(opts,'PERCENT') && opts.PERCENT
    d_str = ['\multicolumn{2}{l}{' Name '} & \multicolumn{2}{c|}{ $\Delta\mathcal{A}_1^e$ (\%)}'];
%     d_str = ['\multicolumn{4}{l|}{' Name '}'];
else
    d_str = ['\multicolumn{4}{l|}{' Name '}'];
end
    [InPeaks_hr,InPeaks_nT] = localMax(Period,InPhase,PeakPeriods);
    [~,OutPeaks_nT] = localMax(Period,Quad,PeakPeriods);
    for ij = 1:length(InPeaks_hr)
        d_str = [d_str ' &' num2str(InPeaks_nT(ij),format)];
        d_str = [d_str ' ' num2str(OutPeaks_nT(ij),format)];
    end
disp([d_str ' \\'])
if opts.INCLUDE_TDs 
    disp('\hline')
end
end %printTableStr
%%
function [out_Hz,out_nT] = localMax(fr,dat,target_fr)
for it = 1:length(target_fr)
    hw = 20; %half-width
    ind = find((abs(fr-target_fr(it))) == min(abs(fr-target_fr(it))));
    inds = ind-hw:ind+hw;
    lmax = find(dat(inds) == max(dat(inds)));
    out_Hz(it) = fr(inds(lmax));
    out_nT(it) = dat(inds(lmax));
end
end %localMax
