function MagResponseParamSpaceJupiter

CALC_NEW = 0;
PLOT_CONTOURS = 1;  % Set to 0 for surfaces
DO_LEGEND = 1;
PLOT_V2020s = 0; 

DO_EUROPA = 1;
DO_GANYMEDE = 1;
DO_CALLISTO = 1;


km = 1e3;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',2*km,'InitialStep',1e-2);

plotText = "";
pos = [0.25, 0.05, 0, 0];

r0 = km;
y0 = 0;
n = 1;

set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultAxesFontSize',15)

npts_k = 50;
npts_D = 60;
np_intp = 200;
    
for ibody = 1:(DO_EUROPA+DO_GANYMEDE+DO_CALLISTO)
    
    if DO_EUROPA && (ibody == 1)
        fprintf("\n @@@ Plotting Europa contours @@@ \n\n")
        
        lname = 'Europa';
        fname = [lname '/' 'EuropaMagParamSpace']; 
        klims = [.1 100];
        Dlims = [1 200];
        
        R_europa_km = 1565;
        R_outer_km = R_europa_km;
        nE = 2*pi/3.55/86400;
        D_ice_km = 20;
        
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
        alevels = { [ 6 9 9.9]
                    [15 50 80 170 190 200 ]
                    [1 4 11 12.5 13.3 14.3 ]};
        plevels = { [10 20 30 40 60 87]
                    [10 30 50 80]
                    [10 20 60 80]};
        
        wsmall = 2*pi/nE./Periods_hr./3600;
        colors = colormap;
        ccolors = colors([10 100 200],:);
        LW = [2 2 2];
        clines = {':','-','-.'};
        
                    color_warmSW = 	[176,0,255]/255;
        V2020.km = [91 117 96 124 91 117 91 119 ];
        V2020.Sm = [0.4132 0.4533 3.3661 3.7646 0.3651 0.3855 2.8862 3.0760];
        V2020.MFCs = {'none','none','b','m','none','none','c',color_warmSW};
        V2020.MECs = {'b','m','k','k','c',color_warmSW,'k','k'};
        V2020.symbols = '^v^v^v^v';
    end
    
    if DO_GANYMEDE && ( (ibody - DO_EUROPA) == 1 )
        fprintf("\n @@@ Plotting Ganymede contours @@@ \n\n")
        
        lname = 'Ganymede';
        fname = [lname '/' 'GanymedeMagParamSpace'];
        klims = [.1 100];
        Dlims = [1 1000];
        
        R_ganymede_km = 2634;
        R_outer_km = R_ganymede_km;
        nE = 2*pi/7.15/86400; % orbit in Hz
        D_ice_km = 50;
        
        % values used in Vance et al. 2020, from Corey Cochrane
        Periods_hr = [171.57 10.53 5.27];
        BnT = [1.21 82.61 2.64];
        alevels = { [0.1 1 1.1]
                    [10 40 70 75 77.5]
                    [1 1.5 2.44 2.47]};
        plevels = { [10 30 60 87]
                    [5 10 30 70 80]
                    [10 20 40 60 80]};
         
        wsmall = 2*pi/nE./Periods_hr./3600;
        colors = colormap;
        ccolors = colors([10 100 200],:);
        LW = [2 2 2];
        clines = {':','-','-.'};
        
        V2020.km = [442 276 458 282];
        V2020.Sm = [0.3890 0.2623 3.1150 1.9483];
        V2020.MFCs = {'none','none','b','m'};
        V2020.MECs = {'b','m','k','k'};
        V2020.symbols = '^v^v';
    end
    
    if DO_CALLISTO && ( (ibody - DO_EUROPA - DO_GANYMEDE) == 1 )
        fprintf("\n @@@ Plotting Callisto contours @@@ \n\n")
        
        lname = 'Callisto';
        fname = [lname '/' 'CallistoMagParamSpace'];
        klims = [.01 100];
        Dlims = [1 1000];
        
        R_callisto_km = 2410.3;
        R_outer_km = R_callisto_km;
        nE = 2*pi/17/86400; % orbit in Hz
        D_ice_km = 50;
        
        % values meant to match Zimmer 2000 
        Periods_hr = [403.04 10.18 3.39];
        BnT = [0.5 50 10]; 
        
        % values used in Vance et al. 2020, from Corey Cochrane
        Periods_hr = [400.33 10.18 5.09];
        BnT = [1.72 37.57 0.25];
        
        alevels = { [0.01 0.05 0.3 0.5 1 1.4 1.5]
                    [12 20 33 35]
                    [0.170 0.226 0.232]};
        plevels = { [ 10 20 40 70]
                    [10 30 60 80]
                    [10 30 60 ]};
                
        wsmall = 2*pi/nE./Periods_hr./3600;
        colors = colormap;
        ccolors = colors([10 100 200],:);
        LW = [2 2 2];
        clines = {':','-','-.'};
        
        V2020.km = [132 21 130 21];
        V2020.Sm = [0.2307 0.0895 1.5256 0.6025];
        V2020.MFCs = {'none','none','b','m'};
        V2020.MECs = {'b','m','k','k'};
        V2020.symbols = '^v^v';

    end

    k_Sm = logspace(log10(klims(1)),log10(klims(2)),npts_k);
    D_ocean_km = logspace(log10(Dlims(1)),log10(Dlims(2)),npts_D);
    
    if CALC_NEW
        Rtop_km = R_outer_km;
        r = linspace(km,Rtop_km*km,5*Rtop_km);

        lk = npts_k;
        lw = length(wsmall);
        ld = npts_D;
        surfnorm = 1;

        [amp,phi,Q]=deal(nan(lk,ld,lw));
        for ik = 1:lk
            parfor io = 1:length(D_ocean_km)
                sig =            [1e-16                             k_Sm(ik)            1e-16        ];
                boundaries = km*[R_outer_km(1)-D_ocean_km(io)-D_ice_km      R_outer_km(1)-D_ice_km R_outer_km(1)  ];
                [~,Q(ik,io,:)] = getMagResponseFunction(r,n,nE,wsmall,sig,boundaries,Rtop_km*km,r0,y0,opts);
                disp([num2str(ik) '/' num2str(lk)  ' of ' num2str(io) '/' num2str(length(D_ocean_km))])
                amp(ik,io,:) = surfnorm*2*abs(Q(ik,io,:));
                phi(ik,io,:) = -angle(Q(ik,io,:))*180/pi;
            end
        end
        save(fullfile([fname 'Dice' num2str(D_ice_km) 'km']));
    else
        sv_opt = [PLOT_CONTOURS, DO_LEGEND, DO_EUROPA, DO_GANYMEDE, DO_CALLISTO];
        ibsv = ibody;
        als = alevels;
        pls = plevels;
        
        
        load(fullfile([fname 'Dice' num2str(D_ice_km) 'km']));
        
        CALC_NEW = 0;
        ibody = ibsv;
        alevels = als;
        plevels = pls;
        PLOT_CONTOURS = sv_opt(1);
        DO_LEGEND = sv_opt(2);
        DO_EUROPA = sv_opt(3);
        DO_GANYMEDE = sv_opt(4);
        DO_CALLISTO = sv_opt(5);
    end

    intMethod = 'makina';

    figno = 111+ibody;
    figure(figno);clf
        set(gcf,'Position',[245   362   852   420]);
    kfine = logspace(log10(k_Sm(1)),log10(k_Sm(end)),np_intp);
    Dfine = logspace(log10(D_ocean_km(1)),log10(D_ocean_km(end)),np_intp);

    subplot(1,2,1);hold on;title({'Amplitude (nT) '})
    set(gca,'xscale','log','yscale','log');box on
    xlabel('$\sigma_\mathrm{ocean}$ (S/m)','Interpreter','latex')
    ylabel('$D_\mathrm{ocean}$ (km)','Interpreter','latex')
    xlim(klims)
    ylim(Dlims)
    
    if PLOT_V2020s
        for ip = 1:length(V2020.km)
            hp = plot(V2020.Sm(ip),V2020.km(ip),V2020.symbols(ip));
            hp.MarkerFaceColor = V2020.MFCs{ip};
            hp.MarkerEdgeColor = V2020.MECs{ip};
        end
    end

    subplot(1,2,2);hold on;title({'Phase Delay (°)'})
    set(gca,'xscale','log','yscale','log');box on
    xlabel('$\sigma_\mathrm{ocean}$ (S/m)','Interpreter','latex')
    %ylabel('$D_\mathrm{ocean}$ (km)','Interpreter','latex')
    xlim(klims)
    ylim(Dlims)
    
    if PLOT_V2020s
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
        ampfine = interp2(k_grid,D_grid,amp(:,:,iw)',k_interp,D_interp,intMethod);
        phifine = interp2(k_grid,D_grid,phi(:,:,iw)',k_interp,D_interp,intMethod);
        if PLOT_CONTOURS
            H(iw) = plotTheContours(kfine,Dfine,BnT(iw),ampfine,phifine,alevels{iw,:},plevels{iw,:},ccolors(iw,:),clines{iw},LW(iw));
        else
            H(iw) = plotSurf(kfine,Dfine,BnT(iw),ampfine,phifine);
            DO_LEGEND = 0;
        end
    end
    
    % Create legend
    if DO_LEGEND
        lstr = strings(1,lw);
        for in = 1:lw
            lstr(in) = ['\begin{tabular}{p{6mm}r}' num2str(Periods_hr(in),'%0.2f') '&hr\end{tabular}'];
        end
        hl = legend(H,lstr,'FontSize',20);
        hl.Interpreter = 'latex';
        hl.Location = 'northeast';
        hl.FontWeight = 'bold';
    end
        
        % Add text label
        ht = text(0.003850892987281,38.68764451136161,0,lname);
        ht.FontSize = 30;
        ht.BackgroundColor = 'w';

    annotation('textbox', pos, 'string', plotText)
    
    fprintf("Printed figure: %d\n", figno)
    
    str2020 = '';
    if PLOT_V2020s
        str2020 = '_WithV2020Models';
    end
    saveas(gcf,['MagPhase' lname str2020 '.eps'],'epsc')
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
end

function SA = plotSurf(k_Sm,D_ocean_km,BnT,amp,phi)
    subplot(1,2,1)
    SA = surf(k_Sm,D_ocean_km,BnT*amp');
    SA.EdgeColor = 'none';
    colorbar
    
    subplot(1,2,2)
    SP = surf(k_Sm,D_ocean_km,phi');
    colorbar
    SP.EdgeColor = 'none';
end