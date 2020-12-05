function MagResponseParamSpaceJupiter

cfg = config;

DO_LEGEND = cfg.DO_LEGEND;

ccolors = [cfg.col_contOrb; cfg.col_contSyn; cfg.col_contHrm];
clines = {cfg.LS_orb, cfg.LS_syn, cfg.LS_hrm};
LW =     [cfg.LW_orb; cfg.LW_syn; cfg.LW_hrm];

km = 1e3;
nbodies = (cfg.DO_EUR + cfg.DO_GAN + cfg.DO_CAL);
r0 = km; y0 = 0; n = 1; opts = cfg.opts_odeParams;
ice_thk.eur = 20;
ice_thk.gan = 50;
ice_thk.cal = 50;

bnames = strings(1,nbodies);
D_ice_km = zeros(1,nbodies);
if cfg.DO_EUR; ei=1; bnames(ei) = 'Europa'; D_ice_km(ei) = ice_thk.eur; end
if cfg.DO_GAN; gi=cfg.DO_EUR+cfg.DO_GAN; bnames(gi) = 'Ganymede'; D_ice_km(gi) = ice_thk.gan; end
if cfg.DO_CAL; ci=cfg.DO_EUR+cfg.DO_GAN+cfg.DO_CAL; bnames(ci) = 'Callisto'; D_ice_km(ci) = ice_thk.cal; end

% Pregenerate figures and get labels
lbl = getPlotLabels(cfg.dft_font, cfg.dft_math, cfg.interpreter);
figs = getMagResFigRefs(lbl, bnames, cfg.NO_PLOTS);

if ~cfg.NO_PLOTS
    for iFig=1:numel(figs.cont)
        set(figs.cont(iFig), 'visible', 'on');
    end
    % Let the user see changes in real time, and only pop up windows at the
    % *start* of the first run on which we generate plots
    drawnow
end

plotText = '';
% Where to place extra text label
pos = [0.25, 0.05, 0, 0];
    
peaks_hr = zeros(nbodies,3);
for ibody = 1:nbodies
    bodyName = char(bnames(ibody));
    fname = [bodyName '/' bodyName 'MagParamSpace'];
    thisSaveStr = [fname 'Dice' num2str(D_ice_km(ibody)) 'km'];
    
    switch char(bodyName)
    case 'Europa'
        fprintf('\n @@@ Plotting Europa contours @@@ \n\n')
        klims = [.1 100];
        Dlims = [1 200];
        fftxlims = [1.5 2e2];

        R_europa_km = 1565;
        R_outer_km = R_europa_km;
        nE = 2*pi/3.55/86400;

        % values to match Khurana 2002
        %BnT = [14 250 20]; 
        %Periods_hr = [85.23 11.23 5.62];
        %alevels = {  [2 11 12.5]
        %            [10 35 210 228 235]
        %            [2 4 8 12 16 18]};
        %plevels = {  [10 20 30 40 60]
        %            []
        %            [20 60 80]};

        % values used in Vance et al. 2020, from Corey Cochrane
        BnT = [10.65 209.78 15.03];
        Periods_hr = [85.20 11.23 5.62];
        peaks_hr(ibody,:) = Periods_hr(:);
        % contour levels
        alevels = { [ 6 9 9.9]
                    [15 50 80 170 190 200 ]
                    [1 4 11 12.5 13.3 14.3 ]};
        plevels = { [10 20 30 40 60 87]
                    [10 30 50 80]
                    [10 20 60 80]};

        V2020.km = [91 117 96 124 91 117 91 119 ];
        V2020.Sm = [0.4132 0.4533 3.3661 3.7646 0.3651 0.3855 2.8862 3.0760];
        V2020.MFCs = {'none','none','b','m','none','none','c',cfg.b000ff};
        V2020.MECs = {'b','m','k','k','c',cfg.b000ff,'k','k'};
        V2020.symbols = '^v^v^v^v';
    
    case 'Ganymede'
        fprintf('\n @@@ Plotting Ganymede contours @@@ \n\n')
        klims = [.1 100];
        Dlims = [1 1000];
        fftxlims = [1 2e3];
        
        R_ganymede_km = 2634;
        R_outer_km = R_ganymede_km;
        nE = 2*pi/7.15/86400; % orbit in Hz
        
        % values used in Vance et al. 2020, from Corey Cochrane
        BnT = [1.21 82.61 2.64];
        Periods_hr = [171.57 10.53 5.27];
        peaks_hr(ibody,:) = Periods_hr(:);
        % contour levels
        alevels = { [0.1 1 1.1]
                    [10 40 70 75 77.5]
                    [1 1.5 2.44 2.47]};
        plevels = { [10 30 60 87]
                    [5 10 30 70 80]
                    [10 20 40 60 80]};
        
        V2020.km = [442 276 458 282];
        V2020.Sm = [0.5166 0.3322 4.0699 2.3476];
        V2020.MFCs = {'none','none',cfg.col_warmestMgSO4,cfg.col_coldestMgSO4};
        V2020.MECs = {cfg.col_warmestMgSO4,cfg.col_coldestMgSO4,'k','k'};
        V2020.symbols = '^v^v';
    
    case 'Callisto'
        fprintf('\n @@@ Plotting Callisto contours @@@ \n\n')
        klims = [.01 100];
        Dlims = [1 1000];
        fftxlims = [1 2e3];
        
        R_callisto_km = 2410.3;
        R_outer_km = R_callisto_km;
        nE = 2*pi/17/86400; % orbit in Hz
        
        % values meant to match Zimmer 2000 
        %Periods_hr = [403.04 10.18 3.39];
        %BnT = [0.5 50 10]; 
        % values used in Vance et al. 2020, from FFTs by Corey Cochrane
        BnT = [1.72 37.57 0.25];
        Periods_hr = [400.33 10.18 5.09];
        peaks_hr(ibody,:) = Periods_hr(:);
        % contour levels
        alevels = {[0.01 0.05 0.3 0.5 1 1.4 1.5]
                   [12 20 33 35]
                   [0.170 0.226 0.232]};
        plevels = {[10 20 40 70]
                   [10 30 60 80]
                   [10 30 60 ]};
                
        V2020.km = [132 21 130 21];
        V2020.Sm = [0.2307 0.0895 1.5256 0.6025];
        V2020.MFCs = {'none','none',cfg.col_warmestMgSO4,cfg.col_coldestMgSO4};
        V2020.MECs = {cfg.col_warmestMgSO4,cfg.col_coldestMgSO4,'k','k'};
        V2020.symbols = '^v^v';
    end

    k_Sm = logspace(log10(klims(1)),log10(klims(2)),cfg.npts_k);
    D_ocean_km = logspace(log10(Dlims(1)),log10(Dlims(2)),cfg.npts_D);
    
    if cfg.CALC_NEW
        r = linspace(km,R_outer_km*km,5*R_outer_km);
        wsmall = 2*pi/nE./Periods_hr./3600;

        lk = cfg.npts_k;
        lw = length(wsmall);
        ld = cfg.npts_D;
        % surfnorm is for normalizing to 1/r^3 dependence, e.g. when
        % R_outer refers to the ionosphere outer boundary. Ae and Q are
        % normalized to max out at the outer conductor boundary, so they
        % must be scaled up if we want to compare to 1.0 being the
        % planetary surface, or scaled down if we were to not model the ice shell
        % in our calculation.
        surfnorm = 1;
        totks = num2str(lk);
        totDs = num2str(ld);

        [amp,phi,Q]=deal(nan(lk,ld,lw));
        % Calculate boundary radii
        io_bdy = R_outer_km - D_ice_km(ibody);
        r_surf = R_outer_km;
        for ik = 1:lk
            parfor io = 1:ld
                % Ocean bottom depth changes as we vary D_ocean
                oc_bot = io_bdy - D_ocean_km(io);
                % Get layered conductivity structure (uniform ocean)
                sig =           [1e-16  k_Sm(ik) 1e-16 ];
                boundaries = km*[oc_bot io_bdy   r_surf];
                % Get Eckhardt (1963) Q calculation
                [~,Q(ik,io,:)] = getMagResponseFunction(r,n,nE,wsmall,sig,boundaries,r_surf*km,r0,y0,opts);
                disp([num2str(ik) '/' totks  ' of ' num2str(io) '/' totDs])
                amp(ik,io,:) = surfnorm*2*abs(Q(ik,io,:));
                phi(ik,io,:) = -angle(Q(ik,io,:))*180/pi;
            end
        end
        save(fullfile(thisSaveStr), 'nbodies', 'k_Sm', 'D_ocean_km', 'amp', 'phi', 'lw');
    else
        load(fullfile(thisSaveStr), 'nbodies', 'k_Sm', 'D_ocean_km', 'amp', 'phi', 'lw');
        if nbodies ~= cfg.DO_EUR + cfg.DO_GAN + cfg.DO_CAL
            error(['ERROR: The file we loaded with CALC_NEW = 0 had the wrong number of bodies saved. ' ...
                'Re-run with CALC_NEW=1. Aborting.'])
        end
    end

    set(0, 'CurrentFigure', figs.cont(ibody));
    clf;
    set(gcf,'Position',[245   362   852   420]);
    kfine = logspace(log10(k_Sm(1)),log10(k_Sm(end)),cfg.np_intp);
    Dfine = logspace(log10(D_ocean_km(1)),log10(D_ocean_km(end)),cfg.np_intp);

    subplot(1,2,1);hold on;title(lbl.ttamp)
    set(gca,'xscale','log','yscale','log');box on
    xlabel(lbl.xcond)
    ylabel(lbl.ythic)
    xlim(klims)
    ylim(Dlims)
    
    if cfg.PLOT_V2020S
        for ip = 1:length(V2020.km)
            hp = plot(V2020.Sm(ip),V2020.km(ip),V2020.symbols(ip));
            hp.MarkerFaceColor = V2020.MFCs{ip};
            hp.MarkerEdgeColor = V2020.MECs{ip};
        end
    end

    subplot(1,2,2);hold on;title(lbl.ttphs)
    set(gca,'xscale','log','yscale','log');box on
    xlabel(lbl.xcond)
    %ylabel(lbl.ythic)
    xlim(klims)
    ylim(Dlims)
    
    if cfg.PLOT_V2020S
        for ip = 1:length(V2020.km)
            hp = plot(V2020.Sm(ip),V2020.km(ip),V2020.symbols(ip));
            hp.MarkerFaceColor = V2020.MFCs{ip};
            hp.MarkerEdgeColor = V2020.MECs{ip};
        end
    end

    H = gobjects(1,lw);
    for iw = 1:lw
        [k_grid,D_grid] = meshgrid(k_Sm,D_ocean_km);
        [k_interp,D_interp] = meshgrid(kfine,Dfine);
        ampfine = interp2(k_grid,D_grid,amp(:,:,iw)',k_interp,D_interp,cfg.intMethod);
        phifine = interp2(k_grid,D_grid,phi(:,:,iw)',k_interp,D_interp,cfg.intMethod);
        if cfg.PLOT_CONTOURS
            H(iw) = plotTheContours(kfine,Dfine,BnT(iw),ampfine,phifine,alevels{iw,:},plevels{iw,:},ccolors(iw,:),clines{iw},LW(iw));
        else
            H(iw) = plotSurf(kfine,Dfine,BnT(iw),ampfine,phifine);
            DO_LEGEND = 0;
        end
    end
    
    % Create legend
    if DO_LEGEND
        lstr = strings(1,lw);
        Hlines = gobjects(1,lw);
        for in = 1:lw
            lstr(in) = ['\begin{tabular}{p{6mm}r}' num2str(Periods_hr(in),'%0.2f') '&hr\end{tabular}'];
            % The following line is a workaround to get Matlab to print
            % straight lines in the legend instead of ellipses for
            % contours.
            Hlines(in) = line([NaN,NaN],[1,1],'Color',ccolors(in,:), 'Linestyle',clines{in}, 'Linewidth',LW(in));
        end
        hl = legend(Hlines,lstr,'FontSize',20);
        hl.Interpreter = 'latex';
        hl.Location = 'northeast';
        hl.FontWeight = 'bold';
    end
        
        % Add text label
        ht = text(0.003850892987281,38.68764451136161,0,bodyName);
        ht.FontSize = 30;
        ht.BackgroundColor = 'w';

    annotation('textbox', pos, 'string', plotText)
    fprintf('Printed contour plot for %s\n', bodyName)
    
    str2020 = '';
    if cfg.PLOT_V2020S
        str2020 = '_WithV2020Models';
    end
    print(figs.cont(ibody),cfg.fig_fmt,fullfile([bodyName '/figures/MagPhase' bodyName str2020 '.eps']));
end

%% get the FFT data processed by Corey Cochrane. 
% Notes:
% all data are in IAU reference frame of reference

% 10 year duration ...
%    data points: 524288
%    Sep 1, 2018 12:00:00 AM UTC
%    Sep 1, 2028 12:00:00 AM UTC
%    Time Resolution: 0.0016611 Hz, 601.9958 seconds, 10.0333 minutes, 0.16722 hours 
%FTdata = ReadGalileanFFTData;
%% set up the peaks to plot and analyze

%% create frequency plots similar to those in Seufert et al. 2011
if cfg.DO_PER
    
    if cfg.DISP_TABLES
        disp('\begin{table}[]')
        disp('\centering')
        disp('    \hline')
        disp('\begin{tabular}{c | c c c}')
        disp('&     \multicolumn{3}{c}{\textbf{Period (hr)}}\\')
        disp('     &     $B_{x,y,z}$ (nT)&     $B_{x,y,z}$ (nT)&     $B_{x,y,z}$ (nT)\\')
        disp('    \hline')
    end
        
    body_name = char(bnames(nbodies));
    FTfile = ['FTdata' body_name '.mat'];
    load(FTfile, 'FTdata');
    Periods =  1./FTdata.frequency./3600; 
    peaks(nbodies) = plotPeriods(figs.ffts(nbodies),body_name,fftxlims,peaks_hr(nbodies,:),FTdata,Periods,lbl);
    for ibody=nbodies-1:-1:1 % Do loop in descending order to avoid preallocating struct
        body_name = char(bnames(ibody));
        FTfile = ['FTdata' body_name '.mat'];
        load(FTfile, 'FTdata');
        peaks(ibody,:) = plotPeriods(figs.ffts(ibody),body_name,fftxlims,peaks_hr(ibody,:),FTdata,Periods,lbl);
    end
    
    if cfg.DISP_TABLES
        for ibody=1:nbodies
            dstr = ['\textbf{' char(bnames(ibody)) '} '];
            for in = 1:3 % count from lowest period (highest frequency)
                dstr = [dstr ' & \textbf{' num2str(peaks(ibody).Y_hr(in),'%0.2f') '}' ];
            end
            disp([dstr ' \\'])
            dstr = '             ';
            for in = 1:3
                dstr = [dstr ' & ' num2str([peaks(ibody).X_nT(in) peaks(ibody).Y_nT(in) peaks(ibody).Z_nT(in)],'%8.2f') ];
            end
            disp([dstr ' \\'])
            disp('    \hline')
        end

        disp('   \end{tabular}')
        disp('    \caption{Caption}')
        disp('     \label{tab:my_label}')
        disp('\end{table}')
    end
end

end

function H = plotTheContours(k_Sm,D_ocean_km,BnT,amp,phi,alevels,plevels,color,line,LW)
    subplot(1,2,1)
    if isempty(alevels)
        [~,H]=contour(k_Sm,D_ocean_km,BnT*amp,'ShowText','on');
    else
        [~,H]=contour(k_Sm,D_ocean_km,BnT*amp,alevels,'ShowText','on');
    end
    H.Color = color;
    H.LineStyle = line;
    H.LineWidth = LW;
    subplot(1,2,2)
    if isempty(plevels)
        [~,H]=contour(k_Sm,D_ocean_km,phi,'ShowText','on');
    else
        [~,H]=contour(k_Sm,D_ocean_km,phi,plevels,'ShowText','on');
    end
    H.Color = color;
    H.LineStyle = line;
    H.LineWidth = LW;
end % plotTheContours

function SA = plotSurf(k_Sm,D_ocean_km,BnT,amp,phi)
    subplot(1,2,1)
    SA = surf(k_Sm,D_ocean_km,BnT*amp');
    SA.EdgeColor = 'none';
    colorbar
    
    subplot(1,2,2)
    SP = surf(k_Sm,D_ocean_km,phi');
    colorbar
    SP.EdgeColor = 'none';
end % plotSurf

function peaks = plotPeriods(fig,body_name,fftxlims,peaks_hr,FTdata,Periods,lbl)
    Bx = FTdata.By_fft; % Swapping x and y FFT amplitudes consistent with
    By = FTdata.Bx_fft; % converting from IAU coordinates to PhiO
    Bz = FTdata.Bz_fft;

    set(0, 'CurrentFigure', fig);
    set(gcf, 'Name', [lbl.fftsp ' ' body_name]);
    clf; hold on;    box on
    set(gcf,'Position',[ 73         207        1097         337])
    set(gca,'xscale','log','yscale','log','xlim',fftxlims,'ylim',[1e-5 300])
    text(1.7,40,body_name,'FontSize',40)
    xlabel(lbl.T_hrs);
    ylabel(lbl.Bcomp)
    plot(Periods,Bx,'b',Periods,By,'k',Periods,Bz,'g');
    
    [peaksX_hr,peaksX_nT] = localMax(Periods,Bx,peaks_hr);
    [peaksY_hr,peaksY_nT] = localMax(Periods,By,peaks_hr);
    [peaksZ_hr,peaksZ_nT] = localMax(Periods,Bz,peaks_hr);
    peaks.X_hr = peaksX_hr;
    peaks.X_nT = peaksX_nT;
    peaks.Y_hr = peaksY_hr;
    peaks.Y_nT = peaksY_nT;
    peaks.Z_hr = peaksZ_hr;
    peaks.Z_nT = peaksZ_nT;
    
    InPhaseX = Bx;
    InPhaseY = By;
    InPhaseZ = Bz;

    QuadX = Bx;
    QuadY = By;
    QuadZ = Bz;

    hp = plot(peaksY_hr,peaksY_nT,'ko');
    hp.MarkerFaceColor = 'k';
    for it = 1:length(peaksY_hr)
        ht = text(peaksY_hr(it),peaksY_nT(it),{[num2str(peaksY_hr(it),'%0.2f') ' hr'],['$B_y$: ' num2str(peaksY_nT(it),'%0.2f') ' nT']});
        ht.Interpreter = 'latex';
        ht.FontSize = 16;
        ht.HorizontalAlignment = 'right';
    end

    hl = legend({[lbl.math 'B_x'],[lbl.math 'B_y'],[lbl.math 'B_z']});
    
end %plotPeriods