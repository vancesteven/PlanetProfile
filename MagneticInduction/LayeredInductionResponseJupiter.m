function outWaveforms = LayeredInductionResponseJupiter(Planet,FTdata,Params)

warning('off','all')
disp('All warnings are turned off. Turn them on again to check for NaN values in the input data.')

cfg = config;

PLOT_WAVEFORMS = cfg.plot_fft;
Planet.PLOT_SIGS = 1;
if ~cfg.calc_new_induc && Planet.PLOT_SIGS; Planet.PLOT_SIGS = 0; end

nTbs = length(Planet.Tb_K);

% Pregenerate figures
lbl = getPlotLabels(cfg.dft_font, cfg.dft_math, cfg.interpreter);
figs = getLayeredFigRefs(lbl, Planet.Tb_K, Planet.name, Planet.PLOT_SIGS, cfg.hold, cfg.no_plots);
set(0,'defaultAxesFontSize',16)

%% increments of distance
km = 1e3;
r0 = km; y0 = 0; n = 1;
opts = cfg.opts_odeLayers;

r = linspace(km,2*Planet.R_m,100);
w = linspace(Params.wlims(1),Params.wlims(2),cfg.npts_w); 
w = exp(w); % cycles per orbit

set(0, 'CurrentFigure', figs.amph);
set(gcf, 'Name', [lbl.amphs ' ' Planet.name]);

if ~cfg.hold; warning('WARNING: cfg.hold=0 is not currently implemented in LayeredInductionResponse. Only one overlayed plot will be printed for each body.'); end
%clf;

for iT = nTbs:-1:1 % Do this loop in descending order to avoid preallocating structs
    fname = [Planet.name 'Waveforms_' char(Planet.Profile_ID(iT))];
    if cfg.calc_new_induc
        disp(['== ' Planet.Ocean.comp ' - ' char(Planet.ice_thk(iT)) ' km'])
        SaveWaveforms = getPlanetMagWaveforms(figs,lbl,r,n,w,Planet.sig(iT,:),Planet.boundaries(iT,:),r0,y0,opts,Planet.f_orb,Planet);
        save(fullfile(Planet.name,fname),'SaveWaveforms')
    else
        try
            load(fullfile(Planet.name,fname),'SaveWaveforms')
        catch loadWaveformsError
            disp(['ERROR: ' fname '.mat was not found.'])
            disp('It probably has not been generated. Set cfg.calc_new_induc to 1 to correct this.')
            rethrow(loadWaveformsError)
        end
    end
    outWaveforms(iT) = SaveWaveforms;
end

%% Figures
if PLOT_WAVEFORMS; plotWaveformData(figs,lbl,outWaveforms,FTdata,Planet,cfg); end

%% Table output
% MJS 2020-10-03:
% Since this portion was written, this function has changed implementation
% dramatically. Previously, a large number of Waveforms files were imported
% and compared, but now we only do 1 composition, 1 body at a time. It
% probably makes sense to cut this, as a working implementation is archived
% with the 1.1.0 release for Vance et al. 2020.
if cfg.disp_tables && cfg.deprecated
    disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
    printWaveformTables(FTdata,outWaveforms,'x');
    disp('YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY')
    printWaveformTables(FTdata,outWaveforms,'y');
    disp('ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ')
    printWaveformTables(FTdata,outWaveforms,'z');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

warning('on','all');
disp('Finished call to LayeredInductionResponse.')
end % LayeredInductionResponse 

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
function outWaveforms = getSw(r,n,w,sig,boundaries,R_m,r0,y0,opts,f_orb,Planet)
    opts.hiprec = true;
    opts.par = false;
%     Q = MagComplexA(boundaries,sig,nE*w,n,opts);
    [k,Q] = getMagResponseFunction(r,n,f_orb,w,sig,boundaries,R_m,r0,y0,opts);
    nbound = length(boundaries(:,1));
    for lm = 1:nbound
        Amplitude(lm,:) = 2*abs(squeeze(Q(:,:,lm)));
        if isfield(Planet,'SURF_NORM') && Planet.SURF_NORM
            Amplitude = Amplitude*(Planet.Rtop/Planet.R_m)^3;
        end
        Phase(lm,:) = -angle(squeeze(Q(:,:,lm)));                
    end
    outWaveforms.Amplitude = Amplitude;
    outWaveforms.Phase = Phase;
    outWaveforms.Q = Q;
    outWaveforms.k = k;
    outWaveforms.w = w;
    outWaveforms.sig = sig;
    outWaveforms.boundaries = boundaries;
end % getSw    
%%
function outWaveforms = getPlanetMagWaveforms(figs,lbl,r,n,w,sig,boundaries,r0,y0,opts,f_orb,Planet)
% get Sw for the different configurations related to the main input model
if Planet.ADD_TRANSITION_BOUNDS
    % Adding discrete transitions makes a more intuitive looking plot of the sigma values. 
    % Setting the transition values to eps seems to change the output InPhase and Quad values by about 0.1%, while setting the value to r0 
    % doesn't seem to change them at all. I am assuming that increments smaller than step size r0 are introducing errors in the ode solver
    btrans = r0; % or eps
    boundaries = [boundaries(1) boundaries(2)-btrans boundaries(2:end-1) boundaries(end-1)+btrans boundaries(end)];  
    sig = [sig(1) eps sig(2:end-1) eps sig(end)];

    b_mionos = [boundaries boundaries(end)+btrans Planet.R_m+Planet.ionos_bounds];
    sPedersen_mionos = [sig Planet.ionosPedersen_sig(1) Planet.ionosPedersen_sig];

    b_ionos = [boundaries boundaries(end)+btrans Planet.R_m+Planet.ionos_bounds];
    sPedersen_ionos = [eps*ones(1,length(sig))  Planet.ionosPedersen_sig(1) Planet.ionosPedersen_sig];
else
    b_mionos = [boundaries Planet.R_m+Planet.ionos_bounds];
    sPedersen_mionos = [sig Planet.ionosPedersen_sig];

    b_ionos = [boundaries Planet.R_m+Planet.ionos_bounds];
    sPedersen_ionos = [1e-16*ones(1,length(sig)) Planet.ionosPedersen_sig];
end

if isfield(Planet,'sub_ocean')
    Dice = Planet.sub_ocean.Depth_km*1e3;
    Dlayer = Planet.sub_ocean.Thickness_km*1e3;
    sig_layer = Planet.sub_ocean.sig;
    b_subo = [boundaries(1)-Dice-Dlayer boundaries(1)-Dice boundaries];
    sPedersen_subo = [1e-16,sig_layer,1e-16,sig(2:end-1),1e-16];
    
    b_subo_ionos = [b_subo Planet.R_m+Planet.ionos_bounds];
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
    
    set(0, 'CurrentFigure', figs.sigs);
    set(gcf, 'Name', [lbl.sigsr ' ' Planet.name])
    clf; box on; hold on;
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
    xlabel(lbl.sigma)
    ylabel(lbl.r_kms)
    set(gca,'xscale','log')
    xlim([1e-6 Inf])
end

Planet.SURF_NORM = false;
disp(' (o)  getting main field');outWaveforms.main = getSw(r,n,w,sig,boundaries,Planet.R_m,r0,y0,opts,f_orb,Planet);

Planet.SURF_NORM = true; 
Planet.Rtop = b_mionos(end);
disp(' +o+  with Pedersen ionsphere');outWaveforms.ionosPedersen = getSw(r,n,w,sPedersen_mionos,b_mionos,b_mionos(end),r0,y0,opts,f_orb,Planet);
if isfield(Planet,'ionosCowling_sig')
    disp('-+o+- with Cowling ionsphere');outWaveforms.ionosCowling = getSw(r,n,w,sCowling_mionos,b_mionos,b_mionos(end),r0,y0,opts,f_orb,Planet);
end
if isfield(Planet,'ionos_only') % only create the ionosphere without the ocean in the few places it's specified among the inputs
    disp(' + +  Pedersen ionsphere only');    outWaveforms.ionos_onlyPedersen = getSw(r,n,w,sPedersen_ionos,b_ionos,b_ionos(end),r0,y0,opts,f_orb,Planet);
    if isfield(Planet,'ionosCowling_sig')
        disp('-+o+- Cowling ionsphere only');outWaveforms.ionos_onlyCowling = getSw(r,n,w,sCowling_ionos,b_mionos,b_mionos(end),r0,y0,opts,f_orb,Planet);
    end
end


Planet.SURF_NORM = false;
disp('  -   mean ocean');outWaveforms.mean = getSw(r,n,w,s_mean,b_mean,Planet.R_m,r0,y0,opts,f_orb,Planet);
disp('  ^   top ocean');outWaveforms.top = getSw(r,n,w,s_top,b_top,Planet.R_m,r0,y0,opts,f_orb,Planet);

if isfield(Planet,'sub_ocean')
    disp(' - -  two oceans');    outWaveforms.sub_ocean = getSw(r,n,w,sPedersen_subo,b_subo,Planet.R_m,r0,y0,opts,f_orb,Planet);
    Planet.SURF_NORM = true; 
    disp(' -+-  two oceans, with ionosphere');    outWaveforms.sub_ocean_ionos = getSw(r,n,w,sPedersen_subo_ionos,b_subo_ionos,b_subo_ionos(end),r0,y0,opts,f_orb,Planet);    
end
end %getPlanetMagWaveforms

%% Output the data as LaTeX tables. The \end{} parts of the tables need to be supplied
function printWaveformTables(FTdata,Waveforms,dir)
makeTableHeader(FTdata,dir)
fnames = fieldnames(Waveforms);
for in = 1:length(fnames)
    plot_opts = getPlotOpts(FTdata.Name,fnames{in},cfg);
    if isfield(plot_opts,'Tb_colorTex')
        table_opts.Tb_colorTex = plot_opts.Tb_colorTex;
    end
    table_opts.dir = dir;
    table_opts.main = Waveforms.(fnames{in}).main; % used for printing PERCENT values
    if isfield(Waveforms.(fnames{in}),'ionos_onlyPedersen')
        table_opts.IONOSPHERE_ONLY = true;
        table_opts.INCLUDE_TDs = false;
        printTableStr('Pedersen',FTdata,Waveforms.(fnames{in}).ionos_onlyPedersen,table_opts);    
        if isfield(Waveforms.(fnames{in}),'ionos_onlyCowling')
            table_opts.IONOSPHERE_ONLY = false;
            printTableStr('Cowling',FTdata,Waveforms.(fnames{in}).ionos_onlyCowling,table_opts);                
        end
        disp(['\hline']) % end the section providing ionosphere strengths
    end
    table_opts.IONOSPHERE_ONLY = false;

    % adiabat--this is the main data source. others can be listed in
    % percent of this value
    table_opts.INCLUDE_TDs = true;
    printTableStr(Waveforms.(fnames{in}).Name,FTdata,Waveforms.(fnames{in}).main,table_opts);
    table_opts.INCLUDE_TDs = false;

    if strcmp(FTdata.Name,'Callisto')
      table_opts.PERCENT = 0;
    else
      table_opts.PERCENT = 1;
    end
        % list the results for the main calculation, adding the ionospheric
    % effects
    printTableStr('Pedersen',FTdata,Waveforms.(fnames{in}).ionosPedersen,table_opts);
    if isfield(Waveforms.(fnames{in}),'ionosCowling')
        printTableStr('Cowling',FTdata,Waveforms.(fnames{in}).ionosCowling,table_opts);                
    end
    
    table_opts.PERCENT = 1;
    %mean ocean
    table_opts.main = Waveforms.(fnames{in}).main;
    printTableStr(['$\overline{\sigma}=$ ' num2str(Waveforms.(fnames{in}).mean.sig(end-3),'%0.4f') ' S/m'],FTdata,Waveforms.(fnames{in}).mean,table_opts);
    %top ocean
    printTableStr(['$\sigma_\mathrm{top}=$ ' num2str(Waveforms.(fnames{in}).top.sig(end-3),'%0.4f') ' S/m'],FTdata,Waveforms.(fnames{in}).top,table_opts);
    table_opts.PERCENT = 0;

    table_opts.PERCENT = 1;
    % any modeled second ocean under high-pressure ice
    if isfield(Waveforms.(fnames{in}),'sub_ocean')
        d_layer = Waveforms.(fnames{in}).sub_ocean.sub_ocean.Thickness_km;
        sig_layer = Waveforms.(fnames{in}).sub_ocean.sub_ocean.sig;
         printTableStr(['\textbf{bottom layer: ' num2str(d_layer) ' km ' num2str(sig_layer) ' S/m}'],FTdata,Waveforms.(fnames{in}).sub_ocean,table_opts); 
         printTableStr('Pedersen',FTdata,Waveforms.(fnames{in}).sub_ocean_ionos,table_opts); 
    end
    table_opts.PERCENT = 0;
    disp('\hline')
end
end %printWaveformTables
%%
function plotWaveformData(figs,lbl,Waveforms,FTdata,Planet,cfg)
for in = 1:length(Waveforms)
    plot_opts = getPlotOpts(Planet.name,Planet.Ocean.comp,Planet.Ocean.w_ocean_pct,Planet.D_Ih_km(in),cfg);
    plotAPhi(figs,lbl,plot_opts,Waveforms(in),FTdata,Planet);
end

% Add labels and save the Aphi figures
savename = [Planet.name '/figures/Ae_' Planet.name];
set(0, 'CurrentFigure', figs.amph);
set(gcf, 'Name', [lbl.amphs ' ' Planet.name]);
print(figs.amph,cfg.fig_fmt,fullfile([savename '_Response' cfg.xtn]));

set(0, 'CurrentFigure', figs.BxAe);
set(gcf, 'Name', ['Bx' lbl.Aecpx ' ' Planet.name]);
labelWaveformFigs('x',Planet.name)
print(figs.BxAe,cfg.fig_fmt,fullfile([savename '_Bx' cfg.xtn]));
set(0, 'CurrentFigure', figs.ByAe);
set(gcf, 'Name', ['By' lbl.Aecpx ' ' Planet.name]);
labelWaveformFigs('y',Planet.name)
print(figs.ByAe,cfg.fig_fmt,fullfile([savename '_By' cfg.xtn]));
set(0, 'CurrentFigure', figs.BzAe);
set(gcf, 'Name', ['Bz' lbl.Aecpx ' ' Planet.name]);
labelWaveformFigs('z',Planet.name)
print(figs.BzAe,cfg.fig_fmt,fullfile([savename '_Bz' cfg.xtn]));

end %plotWaveformData
%% Set the line and marker styles. Conditional formatting strongly depends on the choice and labeling of model inputs. This needs to be generalized for future work.
function plot_opts = getPlotOpts(body_name,o_comp,o_wtpc,ice_thk,cfg)
plot_opts.intp = cfg.intMethod;
plot_opts.fine = cfg.np_wfine;
switch body_name
case 'Europa'
    switch o_comp(1) 
    case 'M'
        plot_opts.line=cfg.ls_Mg;
        if ice_thk > 15
            [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_coldestMgSO4);
            plot_opts.point = '^';
            plot_opts.Tb_colorTex = rgbTex(cfg.col_coldestMgSO4);
        else
            [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_warmestMgSO4);
            plot_opts.point = 'v';
            plot_opts.Tb_colorTex = rgbTex(cfg.col_warmestMgSO4);
        end
        if o_wtpc > 5
            plot_opts.MFC = plot_opts.MEC;
            plot_opts.MEC = 'k';
            plot_opts.LW = cfg.lw_sal;
        else
            plot_opts.MFC = 'none';
            plot_opts.LW = cfg.lw_dil;
        end
    case 'S'
        plot_opts.line=cfg.ls_Sw;
        if ice_thk > 15
            [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_coldestSw);
            plot_opts.point = '^';
            plot_opts.Tb_colorTex = rgbTex(cfg.col_coldestSw);
        else
           [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_warmestSw);
           plot_opts.point = 'v';
           plot_opts.Tb_colorTex = rgbTex(cfg.col_warmestSw);
        end
        if o_wtpc > 5
            plot_opts.MFC = plot_opts.MEC;
            plot_opts.MEC = 'k';
            plot_opts.LW = cfg.lw_sal;
        else
            plot_opts.MFC = 'none';
            plot_opts.LW = cfg.lw_dil;
        end
    end
    % Set zoom level
    plot_opts.Bxxlim = [0 12];
    plot_opts.Bxylim = [0 1.7];
    plot_opts.Byxlim = [0 16];
    plot_opts.Byylim = [0 5];
    plot_opts.Bzxlim = [3 10];
    plot_opts.Bzylim = [0 0.7];
case 'Ganymede'
    plot_opts.line=cfg.ls_Mg;
    if ice_thk > 90
        [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_coldestMgSO4);
        plot_opts.point = '^';
        plot_opts.Tb_colorTex = rgbTex(cfg.col_coldestMgSO4);
    else
        [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_warmestMgSO4);
        plot_opts.point = 'v';
        plot_opts.Tb_colorTex = rgbTex(cfg.col_warmestMgSO4);
    end
    if o_wtpc > 5
        plot_opts.MFC = plot_opts.MEC;
        plot_opts.MEC = 'k';
        plot_opts.LW = cfg.lw_sal;
    else
        plot_opts.MFC = 'none';
        plot_opts.LW = cfg.lw_dil;
    end
    % Set zoom level
    plot_opts.Bxxlim = [0 1.8];
    plot_opts.Bxylim = [0 0.12];
    plot_opts.Byxlim = [0 2.5];
    plot_opts.Byylim = [0 0.6];
    plot_opts.Bzxlim = [0 1.8];
    plot_opts.Bzylim = [0 0.08];
case 'Callisto'
    plot_opts.line=cfg.ls_Mg;
    if ice_thk > 90
        [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_coldestMgSO4);
        plot_opts.point = '^';
        plot_opts.Tb_colorTex = rgbTex(cfg.col_coldestMgSO4);
    else
        [plot_opts.LC,plot_opts.MEC] = deal(cfg.col_warmestMgSO4);
        plot_opts.point = 'v';
        plot_opts.Tb_colorTex = rgbTex(cfg.col_warmestMgSO4);
    end
    if o_wtpc > 5
        plot_opts.MFC = plot_opts.MEC;
        plot_opts.MEC = 'k';
        plot_opts.LW = cfg.lw_sal;
    else
        plot_opts.MFC = 'none';
        plot_opts.LW = cfg.lw_dil;
    end
    % Set zoom level
    plot_opts.Bxxlim = [0 0.008];
    plot_opts.Bxylim = [0 0.01];
    plot_opts.Byxlim = [0 0.2];
    plot_opts.Byylim = [0 0.3];
    plot_opts.Bzxlim = [0 0.002];
    plot_opts.Bzylim = [0 0.012];
end
end %getPlotOpts
%%
function plotDiagnostics(nfig,r,w,PName)
xscale = 'linear';
if ~ishandle(nfig+50)
    figure(nfig+50);
else
    set(0, 'CurrentFigure', nfig+50);
end
subplot(2,2,1);
surf(r/1e3,w,real(k))
set(gca,'xscale',xscale,'zscale','log')
axis tight; 
ylabel(['frequency (cycles per ' PName ' orbit)']);
xlabel('r (km)');
zlabel('Re(k)');

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
function plotAPhi(figs,lbl,plot_opts,Waveforms,FTdata,Planet)
%     global FREQUENCY_AXIS
    f_orb = Planet.f_orb;
    f_small = f_orb/2/pi*Waveforms.main.w; % frequency spectrum
    p_small = 1./f_small./3600;
    
    % Chop off zero point to avoid infinite period
    Period = 1./FTdata.frequency(2:end)/3600;
    Bx = FTdata.By_fft(2:end); % FTdata is in IAU coordinates, so x and y
    By = FTdata.Bx_fft(2:end); % amplitudes are swapped for phiO coordinates.
    Bz = FTdata.Bz_fft(2:end);
    
    LineWidth = plot_opts.LW;

    itype = plot_opts.intp;
    Amp = interp1(p_small,Waveforms.main.Amplitude,Period,itype);
    Phase = interp1(p_small,Waveforms.main.Phase,Period,itype);
    A1e_adiabat = interp1(p_small,Waveforms.main.Q,Period,itype);
    A1e_mean = interp1(p_small,Waveforms.mean.Q,Period,itype);
    A1e_top = interp1(p_small,Waveforms.top.Q,Period,itype);

    ReA1e = Amp.*cos(Phase);
    ImA1e = Amp.*sin(Phase);
%    ReA1e = real(A1e_adiabat);
%    ImA1e = imag(A1e_adiabat);

    InPhaseX = Bx.*ReA1e; 
    InPhaseY = By.*ReA1e;
    InPhaseZ = Bz.*ReA1e;

    QuadX = Bx.*ImA1e; % Same as above, x and y are exchanged here.
    QuadY = By.*ImA1e;
    QuadZ = Bz.*ImA1e;
    
    if ~isfield(plot_opts,'MarkerSize')
        plot_opts.MarkerSize = 10;
    end
    [~,peaksInPhaseX_nT] = localMax(Period,InPhaseX,Planet.peaks_hr);
    [~,peaksInPhaseY_nT] = localMax(Period,InPhaseY,Planet.peaks_hr);
    [~,peaksInPhaseZ_nT] = localMax(Period,InPhaseZ,Planet.peaks_hr);

    [~,peaksQuadX_nT] = localMax(Period,QuadX,Planet.peaks_hr);
    [~,peaksQuadY_nT] = localMax(Period,QuadY,Planet.peaks_hr);
    [~,peaksQuadZ_nT] = localMax(Period,QuadZ,Planet.peaks_hr);
    plot_opts.MS = [4 7 12];
    plot_opts.symbols = '  o';

    plotInductionWaveforms(peaksInPhaseX_nT,peaksQuadX_nT,plot_opts,'x',figs.BxAe);
    plotInductionWaveforms(peaksInPhaseY_nT,peaksQuadY_nT,plot_opts,'y',figs.ByAe);
    plotInductionWaveforms(peaksInPhaseZ_nT,peaksQuadZ_nT,plot_opts,'z',figs.BzAe);

    set(0, 'CurrentFigure', figs.amph);
    subplot(2,1,1);hold on; hl = [];
    hthis = semilogx(Period,Amp,'LineWidth',LineWidth,'LineStyle',plot_opts.line,'Color',plot_opts.LC); 
    xlabel(lbl.T_hrs);
    ylabel(lbl.normA);
    axis tight
    
    ylim([0 1])
    yticks(0:0.2:1)
    hv = vline(Planet.peaks_hr,'r');  
    [hv.Color] = deal('r');
    [hv.LineWidth] = deal(1);
    set(gca,'xscale','log')
    box on

    subplot(2,1,2); hold on
    hpf = semilogx(Period,Phase*180/pi,'LineWidth',LineWidth,'LineStyle',plot_opts.line,'Color',plot_opts.LC);
    xlabel(lbl.T_hrs);
    ylabel(lbl.ttphs)
    axis tight
    ylim([0 90])
    hv = vline(Planet.peaks_hr,'r'); 
    [hv.LineWidth] = deal(1);
    set(gca,'xscale','log')
    box on    
end %plotAPhi
%%
function plotInductionWaveforms(pInPhase,pQuad,plot_opts,dir,nfig)
set(0, 'CurrentFigure', nfig);
set(gcf,'Position',[469   767   589   251])
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
    case 'Ganymede'
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
    case 'Callisto'
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
function printTableStr(Name,FTdata,Planet,opts)
    f_orb = FTdata.f_orb;
    f_small = f_orb/2/pi*Planet.w;
    Period = 1./FTdata.frequency./3600;
    PeakPeriods = 1./FTdata.peaks_Hz/3600;
    
    Amp = interp1(f_small,Planet.Amplitude,FTdata.frequency,'spline');
    Phase = interp1(f_small,Planet.Phase,FTdata.frequency,'spline');

    if isfield(opts,'PERCENT') && opts.PERCENT ==1
        format = '%0.2f'; % output to 0.01% precision
        Amp_main = interp1(f_small,opts.main.Amplitude,FTdata.frequency,'spline');
        Phase_main = interp1(f_small,opts.main.Phase,FTdata.frequency,'spline');
        InPhase = pdiff(Amp.*cos(Phase),Amp_main.*cos(Phase_main));
        Quad = pdiff(Amp.*sin(Phase),Amp_main.*sin(Phase_main));
    else
        format = '%0.3f';
        InPhase = Amp.*cos(Phase).*FTdata.(['B' (opts.dir)]);
        Quad = Amp.*sin(Phase).*FTdata.(['B' (opts.dir)]);  
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
function texCol = rgbTex(rgbVec)
    texCol = ['\color[rgb]{' sprintf("%f,%f,%f",rgbVec) '}'];
end