function outPlanet = PlanetProfile(Planet,Seismic,Params)
%% PlanetProfile
%
% Based on
% S. Vance, M. Bouffard, M. Choukroun, and C. Sotin.
% Ganymede's internal structure including thermodynamics of magnesium sulfate oceans in contact with ice.
% Planetary And Space Science, 96:62-70, 2014.
% (http://dx.doi.org/10.1016/j.pss.2014.03.011)]
%
% Expanded in -- Github Release 1.0
% S. D. Vance, M. P. Panning, S. Staehler, F. Cammarano, B. G. Bills, G. Tobie, S..
% Kamata, S. Kedar, C. Sotin, W. T. Pike, et al.
% Geophysical investigations of habitability in ice-covered ocean worlds.
% Journal of Geophysical Research: Planets, Nov 2018.

%% new function to CheckCompatibility
vernum = PPversion;
disp(['PlanetProfile version ' vernum])
if all(vernum(end-2:end) == 'dev'); disp('This version is in development.'); end
% Check SeaFreeze compatibility and presence on path
% This version of PlanetProfile is compatible with the version number
% below.
seaVer = '0.9.2';
checkSeaFreeze(seaVer);


%% leave as is for now, separate file in python, or nonexistent
% First, get runtime config information
% Fetch this information from an external file so we don't track
% runtime settings in this file
if ~isfield(Params,'cfg')
    cfg = config;
else
    cfg = Params.cfg;
end

%% worry about TauP implementation 
if isfield(Params,'TauP')
    try exist('taupcreate')
        disp ('Taup Model file  will be created. Curve will be plotted. Please use TauP (https://www.seis.sc.edu/taup/) for full functionality')
    catch
        disp('Add mattaup to path')
    end
end
        
       
    %% steve to clean this up by calling Params.cfg
Params.NOPLOTS = cfg.NO_PLOTS;
Params.CALC_NEW = cfg.CALC_NEW;
Params.CALC_NEW_REFPROFILES = cfg.CALC_NEW_REF;
Params.CALC_NEW_SOUNDSPEEDS = cfg.CALC_NEW_SOUND;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY = cfg.CONDUCT;
Params.HOLD = cfg.HOLD;
Params.savefigformat = cfg.fig_fmt;

%% deprecate this!!! steve to remove instances of iT below so additional temperatures will be separate PlanetProfile calls
% Find out how many profiles we will compare
nTbs = length(Planet.Tb_K); 

%% Check for porosity so we know if we will need those plots --> to be moved somewhere more sensible -- make part of CheckFields function, or in python create a class Planet that creates all the fields
POROUS = isfield(Planet,'POROUS_ROCK') && Planet.POROUS_ROCK;

%% CreatePlots
% Show any figures already generated in case user has selected no_plots = 0
% after running with no_plots = 1
if ~Params.NOPLOTS
    figList = findall(0, 'type', 'figure');
    for iFig=1:numel(figList)
        set(figList(iFig), 'visible', 'on');
    end
    % Let the user see changes in real time, and only pop up windows at the
    % *start* of the first run on which we generate plots
    drawnow
end


%% SetupFilenames
POROUS_ICE = isfield(Planet,'POROUS_ICE') && Planet.POROUS_ICE;
if POROUS_ICE; porIceStr = '_PorousIce'; else; porIceStr = ''; end   
% implementing a feature to track silicate composition in output files
% because it's getting confusing as we investigate k2, Q, etc....
if isfield(Seismic,'mantleEOSname'); minEOS = ['_' Seismic.mantleEOSname]; else; minEOS = ''; end
% adjust file name based on keywords(clathrates, porous)
%if isfield(Planet,'Clathrate'); clathStr = '_Clathrates'; else; clathStr = ''; end

savebase = [Planet.name 'Profile_'];
if isfield(Planet,'Clathrate'); savebase = [savebase 'Clathrates_']; end
%savebase = [Planet.name 'Profile_' clathStr];
savefile = [savebase Planet.Ocean.comp ...
    '_' num2str(round(Planet.Ocean.w_ocean_pct)) 'WtPct' minEOS porIceStr ];
    
% MJS 2020-10-02:
% strLow is a placeholder for a calculation that places a relevant string
% into filenames for differentiating mantle compositions. That calculation
% will need to be moved up to account for new additions placing the new
% Profile_fname field into Planet.
strLow = '';
datpath = strcat(Planet.name,'/');
figpath = strcat(Planet.name,'/figures/');

if ~exist(figpath,'dir')
    error(['Figure path ' figpath ' does not exist. Check your Matlab path ' ...
        ' settings and add-with-subfolders the main PlanetProfile directory.'])
end
if ~exist(datpath,'dir')
    error(['Data path ' datpath ' does not exist. Check your Matlab path ' ...
        ' settings and add-with-subfolders the main PlanetProfile directory.'])
end



if ~cfg.SKIP_PROFILES
    
    %% Figure label settings
    set(0,'defaultfigurecolor',[1 1 1]) % white background for figures

    lbl = getPlotLabels(cfg.dft_font, cfg.dft_math, cfg.interpreter);
    % Make shorter versions so they don't muck up the strings too much
    nm = lbl.nm; bnm = lbl.bnm; math = lbl.math;

    %% Figure generation
    % Pre-generate all the figures so that we can work with them in the
    % background, instead of popping up as they are focused.
    % If the windows are still open, we can reuse them
    % without grabbing focus by checking ishandle (or isempty for array values).

    figs = getProfileFigRefs(lbl, Planet.Tb_K, Planet.FeCore, POROUS, cfg.HOLD, cfg.NO_PLOTS);
    % Also grab the layered induction figures in order to pregenerate them
    % since we will be a moment in this function
    if cfg.CALC_NEW_INDUC && ...
            (strcmp(Planet.name,'Europa') || strcmp(Planet.name,'Ganymede') || strcmp(Planet.name,'Callisto'))
        [~] = getLayeredFigRefs(lbl, Planet.Tb_K, Planet.name, Planet.PLOT_SIGS, cfg.HOLD, cfg.NO_PLOTS);
    end
    
    %% Save in a file the densities of adiabats corresponding to different ocean concentrations
    % These are reference profiles, for plotting on hydrosphere density vs
    % pressure plots
    % OceanFreezeDensities TBD
    str_ref_densities = ['ref_densities_' Planet.name '_' Planet.Ocean.comp '_pp.mat'];
    if Params.CALC_NEW_REFPROFILES
        nPr = Params.nsteps_ref_rho;
        %calculate the densities of liquid solution on the liquidus lines
        %corresponding to different concentrations
        wref = Params.wref;
        lw = length(wref);
        Pref_MPa=linspace(0,Params.Pseafloor_MPa,nPr);
        Tref_K = zeros(lw,nPr);
        rho_ref_kgm3 = Tref_K; %allocate
        for jr=1:lw
           for il=1:nPr
              try
                 if ~isfield(Planet.Ocean,'fnTfreeze_K')
                     Tref_K(jr,il) = getTfreeze(Pref_MPa(il),wref(jr),Planet.Ocean.comp);
                 else
                     Tref_K(jr,il) = Planet.Ocean.fnTfreeze_K(Pref_MPa(il),wref(jr));
                 end
              catch
                  disp('caught exception while calculating liquid densities (dashed lines); maybe the search range for fzero is too small?')
                  disp(['il = ' num2str(il)])
                 Tref_K(jr,il) = NaN;
              end
              rho_ref_kgm3(jr,il) = fluidEOS(Pref_MPa(il),Tref_K(jr,il),wref(jr),Planet.Ocean.comp);
           end
           if strcmp(Planet.Ocean.comp,'Seawater')
               isreal = find(~isnan(rho_ref_kgm3(jr,:)));
               rho_ref_kgm3(jr,:)=interp1(Pref_MPa(isreal),rho_ref_kgm3(jr,isreal),Pref_MPa,'linear','extrap');
           end
        end
        save(fullfile([datpath str_ref_densities]),'rho_ref_kgm3','Pref_MPa','Tref_K')
    else
        try
            load(fullfile([datpath str_ref_densities]));
        catch
            error(['ERROR: A reference density file for ' Planet.name ' ' Planet.Ocean.comp ' was not found. Re-run with calc_new_ref=1 to generate the needed file.'])
        end
    end
    
end % ~cfg.SKIP_PROFILES

%% globals for functions in fzero --> just swEOS_chooser this needs to be cleaned up --> SetupEOS
    wo = Planet.Ocean.w_ocean_pct;
if strcmp(Planet.Ocean.comp,'Seawater')
    global swEOS
    swEOS.gsw = swEOS_chooser('gsw305');
elseif strcmp(Planet.Ocean.comp,'NaCl')
    error(['NaCl is not currently implemented.'])
elseif strcmp(Planet.Ocean.comp,'NH3')        
%     error(['NH3 is not currently implemented.'])
elseif strcmp(Planet.Ocean.comp,'MgSO4')
    conduct_scaling_MgSO4 = (1+4*wo); % empirical scaling of electrical conductivity from 1 bar values compiled in Hand and Chyba 2007
end

if isfield(Seismic,'mantleEOS')
    thiseos = split(Seismic.mantleEOS,'.tab');
    thiseos = char(thiseos(1));
else
    thiseos = 'none';
end
bar2GPa = 1e-4;

Gg = 6.67300e-11; % m3 kg-1 s-2


% Filename strings
vsP = 'Porosity_vs_P';
vsR = 'Porosity_vs_R';
vperm = 'Permeability';
vgsks = 'Gs_Ks';
vseis = 'Seismic';
vcond = 'Conductivity';
vgrav = 'Gravity';
vmant = 'MantleDens';
vcore = 'CoreMantTrade';
vpvt4 = 'PTx4';
vpvt6 = 'PTx6';
vwedg = 'Wedge';


%% FEATURE!  COMING SOON!
if ~isfield(Planet,'NoH2O') % backward compatibility--haven't finished implementing water-free worlds
    Planet.NoH2O =0;
end

%SetupClathrates
[n_clath, n_iceI, n_ocean] = deal(zeros(1,nTbs));
if isfield(Planet,'Clathrate')
disp('Running with clathrate parameters')
if isfield(Params,'Clathrate_set_depth') % checks if clathrates have set depth
    max_clath_depth=Params.Clathrate_set_depth;
    Check_clath_depth=1;
else
    max_clath_depth=10e15;
    Check_clath_depth=0;
end
else
Params.nsteps_clath = 0;
Check_clath_depth=0;
max_clath_depth=1e15;% ridiculous high number
% MJS 2021-10-27: This format is being removed in the python
% implementation. Using the already-included logical flag
% Planet.CLATHRATE as a check is a much more sensible way to
% toggle clathrate modeling.
end
%% SetupLayers
nsteps = Params.nsteps_iceI + Params.nsteps_ocean + Params.nsteps_clath;

% MJS 2021-10-28: This got moved up from below where saving/reloading a la
% CALC_NEW was done. We needed to do a bit more calculations before making
% the save files.
[z_m,r_m,g_ms2,M_above_kg,M_below_kg] = deal(zeros(nTbs,nsteps));
M_above_kg(:,1) = 0;
M_below_kg(:,1) = Planet.M_kg;
r_m(:,1) = Planet.R_m;
g_ms2(1:nTbs,1) = Gg*Planet.M_kg/Planet.R_m^2;
D_conductivityIh = 632; % W m-1; Andersson et al. 2005 (For comparison, Mckinnon 2006 uses a value of 621 from Slack 1980)

[Planet.Profile_fname, Planet.Profile_ID] = deal(strings(1,nTbs));

% Preallocate for ConvectionDeschampsSotin
[Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I,Qb,Planet.zb_outerIce_m] = deal(zeros(1,nTbs));


%%
if ~Planet.NoH2O
    if Params.CALC_NEW
        % adds number of clathrates based on inputs or sets to zero, also sets
        % maximum depth if specified

%%      

        [T_K,P_MPa,rho_kgm3] = deal(zeros(nTbs,nsteps));

        phase = zeros(nTbs,nsteps);
        % assign phase based on clathrates or Ih
        phase(:,(1+Params.nsteps_clath):(Params.nsteps_clath+Params.nsteps_iceI))=1;
        phase(:,1:Params.nsteps_clath)=30; % index to indicate clathrates, which are not an ice phase
        Tb_K = Planet.Tb_K;

        %preallocate. Required to test for max clathrate depth
        z_m = zeros(nTbs,nsteps);
        g_ms2 = z_m;
        [M_above_kg,M_below_kg] = deal(rho_kgm3); % mass above and below the silicate interface
        M_above_kg(:,1) = 0;
        M_below_kg(:,1) = Planet.M_kg;
        [z_m,r_m] = deal(zeros(nTbs,nsteps));
        r_m(:,1) = Planet.R_m;
        g_ms2(1:nTbs,1) = Gg*Planet.M_kg/Planet.R_m^2;
    

        %% IceLayers --- function iteratively setting up the thermal profile, the density and temperature of the layer with each pressure step
        %ice Ih, ice III, ice V
        %--------------------------------------------------------------------------
        for iT = 1:nTbs %draw thermal profiles corresponding to the different choices of temperature at the bottom of the Ice I shell
            
            T_K(iT,1) = Planet.Tsurf_K;
            P_MPa(iT,1) = Planet.Psurf_MPa;
            % if ice shell had to be thinned in previous run, this will
            % reset indexing correctly. Set starting values:
            n_iceI(iT) = Params.nsteps_iceI; % these are redundant. just change text below to use Params.n...
            n_clath(iT) = Params.nsteps_clath;
            n_ocean(iT) = Params.nsteps_ocean; % 

            if n_clath(iT)>0
                clath_out{iT,1} = Helgerud_sI(P_MPa(iT,1),T_K(iT,1));
                rho_kgm3(iT,1)=clath_out{iT,1}.rho;
                [Cp(iT,1) alpha_K(iT,1)]= getCpIce(P_MPa(iT,1),T_K(iT,1),phase(iT,1)) ;
            else
                rho_kgm3(iT,1) = getRhoIce(P_MPa(iT,1),T_K(iT,1),1);
                [Cp(iT,1) alpha_K(iT,1)]= getCpIce(P_MPa(iT,1),T_K(iT,1),phase(iT,1)) ;
            end
            try
                Pb_MPa(iT) = getPfreeze(Tb_K(iT),wo,Planet.Ocean.comp);
                %[Cp(iT,1) alpha_K(iT,1)]= getCpIce(P_MPa(iT,1),T_K(iT,1),phase(iT,1)) ;


                deltaP = Pb_MPa(iT)/(n_clath(iT)+n_iceI(iT)-1);
                % assumes Pressure gradient doesnt change betweeen clathrates
                % and need to check assumption
                % if n_clath == 0 this will skip and move on to next section
                %% ClathrateLayer function within IceLayers
                for il=2:n_clath(iT) % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                    P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                    T_K(iT,il) = (Planet.Tb_K(iT).^(P_MPa(iT,il)./Pb_MPa(iT))).*(Planet.Tsurf_K.^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                    clath_out{iT,il} = Helgerud_sI(P_MPa(iT,il),T_K(iT,il));

                    rho_kgm3(iT,il)=clath_out{iT,il}.rho;
                    [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,il)) ;

                    z_m(iT,il) = z_m(iT,il-1)+ (P_MPa(iT,il)-P_MPa(iT,il-1))*1e6/g_ms2(iT,il-1)/rho_kgm3(iT,il-1);
                    r_m(iT,il) = Planet.R_m-z_m(iT,il);

                    % determine local gravity and check to make sure maximum
                    % clathrate depth isn't reached
                    M_below_kg(iT,il) = M_above_kg(iT,il-1) + 4/3*pi*(r_m(iT,il-1)^3-r_m(iT,il)^3)*rho_kgm3(iT,il);
                    M_below_kg(iT,il) = Planet.M_kg-M_above_kg(iT,il);
                    g_ms2(iT,il) = Gg*M_below_kg(iT,il)/r_m(iT,il)^2;
                    if z_m(iT,il)>max_clath_depth % check to see if maximum depth of clathrates reached
                        n_iceI(iT)=n_iceI(iT)+(n_clath(iT)-il)+1;
                        n_clath(iT)=il-1;
                        phase(:,(1+n_clath(iT)):(n_clath(iT)+n_iceI(iT)))=1;
                        break
                    end
                    %rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
                end % end ClathrateLayer
                if n_clath(iT)>0
                   ii=n_clath(iT)+1;
                else
                    ii=2;
                end
                
                %% IceILayer
                for il=ii:(n_iceI(iT)+n_clath(iT)) % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                    P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                    T_K(iT,il) = (Planet.Tb_K(iT).^(P_MPa(iT,il)./Pb_MPa(iT))).*(Planet.Tsurf_K.^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                    rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1); % I THINK THIS CALL CAN ACCEPT A VECTOR
                    [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,il)) ; % I THINK THIS CALL CAN ACCEPT A VECTOR

                    %rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
                end


                nIceIIILithosphere=0; 
                nIceVLithosphere=0; 
                PbI_MPa = Pb_MPa;
            catch
                disp('PlanetProfile failed to get Pb! Here''s why:')
                disp(lasterr)
                disp('Maybe it''s okay? If execution stopped, then probably not. Try looking where getPfreeze is called.')

                %% IceIIIUnderplateLayer
                if isfield(Params,'BOTTOM_ICEIII') && Params.BOTTOM_ICEIII % this will elicit an error if one has set the temperature too low but one hasn't specified ice III or V under the ice I
                    disp('Adding ice III to the bottom of ice Ih. Make sure the ocean salinity is high enough that doing this makes sense')
                    nIceIIILithosphere=5;
                    phase(iT,(n_iceI(iT)+n_clath(iT))-nIceIIILithosphere:(n_iceI(iT)+n_clath(iT)))=3;

                    PbI_MPa(iT) = 210; % the Ih-III transition is essentially fixed, within a few MPa
                    deltaP = PbI_MPa(iT)/((n_iceI(iT)+n_clath(iT))-5-1);  % save five entries at the bottom for ice III

                    for il=2:n_clath(iT) % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                    P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                    T_K(iT,il) = (Planet.Tb_K(iT).^(P_MPa(iT,il)./Pb_MPa(iT))).*(Planet.Tsurf_K.^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                    clath_out{iT,il} = Helgerud_sI(P_MPa(iT,il),T_K(iT,il));

                    rho_kgm3(iT,il)=clath_out{iT,il}.rho;
                     z_m(iT,il) = z_m(iT,il-1)+ (P_MPa(iT,il)-P_MPa(iT,il-1))*1e6/g_ms2(iT,il-1)/rho_kgm3(iT,il-1);
                     r_m(iT,il) = Planet.R_m-z_m(iT,il);
                     [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,il))

                     % determine local gravity
                     M_above_kg(iT,il) = M_above_kg(iT,il-1) + 4/3*pi*(r_m(iT,il-1)^3-r_m(iT,il)^3)*rho_kgm3(iT,il);
                     M_below_kg(iT,il) = Planet.M_kg-M_above_kg(iT,il);
                     g_ms2(iT,il) = Gg*M_below_kg(iT,il)/r_m(iT,il)^2;
                     if z_m(iT,il)>max_clath_depth % check to see if maximum depth of clathrates reached
                        n_iceI(iT)=n_iceI(iT)+(n_clath(iT)-il)+1;
                        n_clath(iT)=il-1;
                        phase(:,(1+n_clath(iT)):(n_clath(iT)+n_iceI(iT)))=1;
                        break
                     end 

                    %rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
                    end
                    if n_clath(iT)>0
                        ii=n_clath(iT)+1;
                    else
                        ii=2;
                    end

                    for il=ii:(n_iceI(iT)+n_clath(iT))-nIceIIILithosphere % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                       P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                     T_K(iT,il) = (Planet.Tb_K(iT).^(P_MPa(iT,il)./Pb_MPa(iT))).*(Planet.Tsurf_K.^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                        rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
                        [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,il)) ;

                    %rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
                    end


                    TbIII = Tb_K(iT)+1;
                    Pb_MPa(iT) = getPfreezeIII(TbIII,wo,Planet.Ocean.comp);% fix the thickness of the ice III layer based on T>Tb by 2K
                    deltaP = (Pb_MPa(iT)-PbI_MPanIceIIILithosphere);  %  five entries at the bottom for ice III
                    for il = (n_clath(iT)+n_iceI(iT))-nIceIIILithosphere+1:(n_iceI(iT)+n_clath(iT))
                        P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                        T_K(iT,il) = (TbIII.^(P_MPa(iT,il)./Pb_MPa(iT))).*(TbIII(iT).^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                        rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),3);
                        [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),3) ;
                    end

                elseif isfield(Params,'BOTTOM_ICEV') && Params.BOTTOM_ICEV
                    disp('Adding ice V and III to the bottom of ice Ih. Make sure the ocean salinity is high enough that doing this makes sense')
                    nIceIIILithosphere=5;
                    nIceVLithosphere=5;
                    phase(iT,n_iceI(iT)-nIceVLithosphere-nIceIIILithosphere:n_iceI(iT)-nIceVLithosphere)=3;
                    phase(iT,n_iceI(iT)-nIceVLithosphere:n_iceI(iT))=5;

                    PbI_MPa(iT) = 210; % the Ih-V transition is essentially fixed, within a few MPa
                    deltaP = PbI_MPa(iT)/((n_iceI(iT)+n_clath(iT))-5-1);  % save five entries at the bottom for ice V

                    %clathrates
                    for il=2:n_clath(iT) % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                          P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                          T_K(iT,il) = (Planet.Tb_K(iT).^(P_MPa(iT,il)./Pb_MPa(iT))).*(Planet.Tsurf_K.^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                          clath_out{iT,il} = Helgerud_sI(P_MPa(iT,il),T_K(iT,il));

                          rho_kgm3(iT,il)=clath_out{iT,il}.rho;
                          [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,il))
                          z_m(iT,il) = z_m(iT,il-1)+ (P_MPa(iT,il)-P_MPa(iT,il-1))*1e6/g_ms2(iT,il-1)/rho_kgm3(iT,il-1);
                          r_m(iT,il) = Planet.R_m-z_m(iT,il);

                          % determine local gravity
                          M_above_kg(iT,il) = M_above_kg(iT,il-1) + 4/3*pi*(r_m(iT,il-1)^3-r_m(iT,il)^3)*rho_kgm3(iT,il);
                          M_below_kg(iT,il) = Planet.M_kg-M_above_kg(iT,il);
                          g_ms2(iT,il) = Gg*M_below_kg(iT,il)/r_m(iT,il)^2;
                          if z_m(iT,il)>max_clath_depth % check to see if maximum depth of clathrates reached
                              n_iceI(iT)=n_iceI(iT)+(n_clath(iT)-il)+1;
                              n_clath(iT)=il-1;
                              phase(:,(1+n_clath(iT)):(n_clath(iT)+n_iceI(iT)))=1;
                              break
                        end

                        %rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
                    end
                    if n_clath(iT)>0
                        ii=n_clath(iT)+1;
                    else
                        ii=2;
                    end
                    
                    %ice Ih
                    % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                    for il=ii:(n_iceI(iT)+n_clath(iT))-nIceVLithosphere % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                        P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                        T_K(iT,il) = (Planet.Tb_K(iT).^(P_MPa(iT,il)./Pb_MPa(iT))).*(Planet.Tsurf_K.^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                        rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),phase(iT,ill));
                        [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,ill)) ;
                        %rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
                    end
                    
                    %ice III
                    TbIII = Tb_K(iT)+1;
                    Pb_MPa(iT) = getPfreezeIII(TbIII,wo,Planet.Ocean.comp);% fix the thickness of the ice III layer based on T>Tb by 2K
                    deltaP = (Pb_MPa(iT)-PbI_MPa)/(nIceIIILithosphere);  %  five entries at the bottom for ice III
                    for il = (n_iceI(iT)+n_clath(iT))-nIceIIILithosphere+1:(n_iceI(iT)+n_clath(iT))
                        P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                        T_K(iT,il) = (TbIII.^(P_MPa(iT,il)./Pb_MPa(iT))).*(TbIII(iT).^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                        rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),3);
                        [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),3) ;
                    end

                    %ice V
                    TbV = Tb_K(iT)+1;
                    Pb_MPa(iT) = getPfreezeV(TbV,wo,Planet.Ocean.comp);% fix the thickness of the ice V layer based on T>Tb by 2K
                    deltaP = (Pb_MPa(iT)-PbI_MPa)/(nIceVLithosphere);  %  five entries at the bottom for ice V
                    for il = (n_iceI(iT)+n_clath(iT))-nIceVLithosphere+1:(n_iceI(iT)+n_clath(iT))
                        P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
                        T_K(iT,il) = (TbV.^(P_MPa(iT,il)./Pb_MPa(iT))).*(TbV(iT).^(1-P_MPa(iT,il)./Pb_MPa(iT)));
                        rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),3);
                        [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),3) ;
                    end

                    %adjusts for surface porosity if set
                    if POROUS_ICE % Steve to ask Angela if this is actually being used.
                        % correction for porosity
                        por_in(iT).p=P_MPa(iT,1:il)*1e-3;
                        por_in(iT).t = T_K(iT,1:il);
                        por_in(iT).den = rho_kgm3(iT,1:il);
                        %por_in(iT).vp = velsIce.Vclathl_kms(iT,:);
                        %por_in(iT).vs = velsIce.Vclatht_kms(iT,:);
                        if isfield(Planet,'phi_surface')
                            por_out(iT) = get_porosity_ice(por_in(iT),Planet.phi_surface);
                        else
                            por_out(iT) =get_porosity_ice(por_in(iT));
                        end
                        %permeability = por_out(iT).per;
                        %ice_ind=find(phase(iT,:)==30);
                        rho_kgm3(iT,1:il) = por_out(iT).den;
                        por(iT,1:il)=por_out(iT).por;
                        %velsIce.Vclathl_kms(iT,:) = por_out(iT).vp;
                        %velsIce.Vclatht_kms(iT,:) = por_out(iT).vs;

                        disp(['Average Porosity: ' num2str(mean(por_out(iT).por))])
                        disp(['Porosity: ' num2str(por_out(iT).por)])
                    end
                end % end IceIIIUnderplateLayer
            end

                %% OceanLayer
            %OCEAN + ICE III/V/VI SHELL
            %--------------------------------------------------------------------------
            deltaP = (Params.Pseafloor_MPa-Pb_MPa(iT))/n_ocean(iT); %
            if deltaP<=0 % prevent an embarassing error that can occur when adapting a file for a small object to that of a larger one.
                error('negative increment of pressure while computing ocean profile. Be sure Params.Pseafloor_MPa is greater than the likely pressure at the silicate interface. Currently it''s less than the pressure at the base of the ice I layer.')
            end
            for il =1:n_ocean(iT)
                ill = il+n_iceI(iT)+n_clath(iT);
                disp(['iT: ' num2str(iT) '; il: ' num2str(il) '; P_MPa: ' num2str(round(P_MPa(iT,ill-1))) '; T_K: ' num2str(round(T_K(iT,ill-1)))]);
                P_MPa(iT,ill) = Pb_MPa(iT) + il*deltaP;

                if il==1 % establish the phase vector
                    phase(iT,ill) = getIcePhase(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
                elseif phase(iT,ill-1)~=6
                    phase(iT,ill) = getIcePhase(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
                    if phase(iT,ill-1)==3 && phase(iT,ill)==0 % fix to instabilities in the phase lookup for ice III
                        disp('Fixing apparent instability in the EOS?ice III is forming above the ocean where it wasn''t specified')
                        phase(iT,ill)=3;
                    end
                    if phase(iT,ill-1)==0 && phase(iT,ill)==1 % fix to instabilities in the phase lookup for ice I
                        disp('Fixing apparent instability in the EOS?ice Ih is forming under the ocean')
                        phase(iT,ill)=0;
                    end
                    if phase(iT,ill-1)==5 && phase(iT,ill)==0 % fix to instabilities in the phase lookup for ice V for Titan NH3 LBF svance Jul 6 2021
                        disp('Fixing apparent instability in the EOS? fluid is forming under ice V')
                        phase(iT,ill)=5;
                    end
                else
                    phase(iT,ill) = 6;
                end

                if phase(iT,ill) == 0 %if ocean
                    [rho_ocean,Cp(iT,ill),alpha_o]= fluidEOS(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);
                    % Forbidding MgSO4 in the check below avoids a problem in
                    % application of the EOS for low-salinity MgSO4 oceans.
                    % Negative thermal expansivity regions in the MgSO4 EOS may be
                    % artifacts of the current EOS calculation. See Vance et al. 2014
                    if alpha_o<=0 && ~strcmp(Planet.Ocean.comp,'MgSO4') && ~Planet.ALLOW_NEGALPHA
                        disp('Ocean alpha at ice interface is less than zero. adjusting temperature upward.')
                        disp('This means there''s a conductive layer at the interface with thickness inversely proportional to the heat flow.')
                        disp('The thickness is likely less than a few 100 m. See Melosh et al. 2004.')
                        Tnew = fzero(@(T) alphaAdjust(P_MPa(iT,ill),T,wo,Planet.Ocean.comp),273);
                        disp(['The new temperature is ' num2str(Tnew) ', or delta T = ' num2str(Tnew-T_K(iT,ill-1)) '.'])
                        T_K(iT,ill)=Tnew;
                    else
                        T_K(iT,ill) = T_K(iT,ill-1)+ alpha_o*T_K(iT,ill-1)./(Cp(iT,ill))./(rho_ocean)*deltaP*1e6; %adiabatic gradient in ocean; this introduces an error, in that we are using the temperature from the previous step
                    end
                    rho_kgm3(iT,ill) = fluidEOS(P_MPa(iT,ill),T_K(iT,ill),wo,Planet.Ocean.comp); % to be deleted
                    
%                         [rho_ocean,Cp(iT,ill),alpha_o]= fluidEOS(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);% THIS IS REDUNDANT TO THE CALCULATION OF RHO_OCEAN ABOVE
                else %% HPIceLayers
                    % allow for stable dense fluid under high pressure ice --> this
                    % was added for the callisto study. Commented out currently
                    rhoice = getRhoIce(P_MPa(iT,ill),T_K(iT,ill-1),phase(iT,ill));
                    % rhoocean = fluidEOS(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);
                    % if rhoocean>=rhoice
                    %     phase(iT,ill)=0;
                    % end
                    if ~isfield(Planet.Ocean,'fnTfreeze_K')
                        if il==6
                            x = 1;
                        end
                        T_K(iT,ill) = getTfreeze(P_MPa(iT,ill),wo,Planet.Ocean.comp,T_K(iT,ill-1)); %find the temperature on the liquidus line corresponding to P(k,i+nIceI); should really use a conductive profile here, but that would seem to naturally bring us back to the liquidus. Steve wants to confirm.
                    else
                        T_K(iT,ill) = Planet.Ocean.fnTfreeze_K(P_MPa(iT,ill),wo);
                    end
                    if T_K(iT,ill)<T_K(iT,ill-1)
                        T_K(iT,ill)=T_K(iT,ill-1); % this may no longer be needed because negalpha is accounted for above -- sincerely, steve 6/15/21
                    end
                    rho_kgm3(iT,ill) = getRhoIce(P_MPa(iT,ill),T_K(iT,ill),phase(iT,ill));

%                     [rho_ocean,Cp(iT,ill),alpha_o]=
%                     fluidEOS(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);
%                     % this is odd and doesn't seem to be needed--steve
%                     6/15/21
                    %[Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,ill)) ;
                end
            end
            
            % MJS 2021-10-29: I don't think these lines are important or do
            % anything substantial. Both quantities are recalculated
            % aplenty in the next loops. This was being done just before
            % saving/reloading in the master branch before now.
            rho_kgm3(iT,1) = rho_kgm3(iT,2); % continuity
            Cp(iT,1)=Cp(iT,2);
    

    %%%%%%%%%%%%%%%%%%%%%
    % convert to depth —- PlanetDepths
    %%%%%%%%%%%%%%%%%%%
    %% calculate gravity in each layer instead of assuming surface gravity applies.
    % allocate variables
    %% HydrosphereDepths
        deltaP = Pb_MPa(iT)/(n_iceI(iT)+n_clath(iT));

        % calculates depth for clathrates and ice separatly so the depths
        % can be recorded accurately
        for il = 2:(n_clath(iT))
            % calculate depth
            z_m(iT,il) = z_m(iT,il-1)+ (P_MPa(iT,il)-P_MPa(iT,il-1))*1e6/g_ms2(iT,il-1)/rho_kgm3(iT,il-1);
            % convert to radius
            r_m(iT,il) = Planet.R_m-z_m(iT,il); 

            % determine local gravity
            M_above_kg(iT,il) = M_above_kg(iT,il-1) + 4/3*pi*(r_m(iT,il-1)^3-r_m(iT,il)^3)*rho_kgm3(iT,il);
            M_below_kg(iT,il) = Planet.M_kg-M_above_kg(iT,il);
            g_ms2(iT,il) = Gg*M_below_kg(iT,il)/r_m(iT,il)^2;
        end
        %checks if there were clathrates or not
        if n_clath(iT)>0
            ii=n_clath(iT)+1;
            Zclath(iT)=z_m(iT,n_clath(iT));
            Planet.Zclath(iT)=Zclath(iT);
        else
            ii=2;
            Zclath(iT)=0;
            Planet.Zclath(iT)=0;
        end
        for il = ii:(n_clath(iT)+n_iceI(iT))
            % calculate depth
            z_m(iT,il) = z_m(iT,il-1)+ (P_MPa(iT,il)-P_MPa(iT,il-1))*1e6/g_ms2(iT,il-1)/rho_kgm3(iT,il-1);
            % convert to radius
            r_m(iT,il) = Planet.R_m-z_m(iT,il); 

            % determine local gravity
            M_above_kg(iT,il) = M_above_kg(iT,il-1) + 4/3*pi*(r_m(iT,il-1)^3-r_m(iT,il)^3)*rho_kgm3(iT,il);
            M_below_kg(iT,il) = Planet.M_kg-M_above_kg(iT,il);
            g_ms2(iT,il) = Gg*M_below_kg(iT,il)/r_m(iT,il)^2;
        end
        if isempty(il)
            zb_outerIce_m(iT)=z_m(iT,ii-1);
            Planet.zb_outerIce_m(iT)=zb_outerIce_m(iT);
        else
            zb_outerIce_m(iT)=z_m(iT,il);
            Planet.zb_outerIce_m(iT)=zb_outerIce_m(iT);
        end
        deltaP = (Params.Pseafloor_MPa-Pb_MPa(iT))/n_ocean(iT); %
        for il = 1:n_ocean(iT)
            ill = il+(n_iceI(iT)+n_clath(iT));
            %calculate depth
            dz = deltaP*1e6/g_ms2(iT,il-1+(n_iceI(iT)+n_clath(iT)))/rho_kgm3(iT,ill);
            z_m(iT,ill) = z_m(iT,il-1+(n_iceI(iT)+n_clath(iT)))+ dz; % using the previous gravity step, since we haven't calculated the present step.  this introduces an error
            % convert to radius
            r_m(iT,ill) = Planet.R_m-z_m(iT,ill); 

            % determine local gravity
            M_above_kg(iT,ill) = M_above_kg(iT,il-1+(n_iceI(iT)+n_clath(iT)))+4/3*pi*((r_m(iT,ill)+dz)^3-r_m(iT,ill)^3)*rho_kgm3(iT,ill);
            M_below_kg(iT,ill) = Planet.M_kg-M_above_kg(iT,ill);
            g_ms2(iT,ill) = Gg*M_below_kg(iT,ill)/r_m(iT,ill)^2;
        end
        disp(['z_iceI: ' num2str(zb_outerIce_m(iT)/1e3) ' km'])
        z_ocean_m(iT)= z_m(iT,n_ocean(iT)); % depth to the ocean
        
        

        %% compute conductive heat through the ice I layer
        Qb(iT) = D_conductivityIh*log(Planet.Tb_K(iT)/Planet.Tsurf_K)/Planet.zb_outerIce_m(iT);
        
        %% compute solid state convection ice
        % We use these values much later, but we need them to construct
        % the filename shortly, so we calculate them now. This allows
        % us to skip unnecessary calculations in the case of
        % cfg.SKIP_PROFILES=1
          if max_clath_depth<Planet.zb_outerIce_m(iT) % only a clathrate lid
            % asummes Q across ice-ocean is same Q across clathrates/ice. 
           [Q_Wm2(iT), T_clath_ice(iT), deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),nu(iT),CONVECTION_FLAG_I(iT)]= ...
       clathrate_lid_thermo(Planet.Tsurf_K,Planet.Tb_K(iT),P_MPa(iT,1:n_clath(iT)+n_iceI(iT)), n_clath(iT),n_iceI(iT),Planet.zb_outerIce_m(iT), max_clath_depth,g_ms2(iT,1));      
          else
           [Q_Wm2(iT),deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),nu(iT),CONVECTION_FLAG_I(iT)]=...
            ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K(iT),PbI_MPa(iT)/2,Planet.zb_outerIce_m(iT),g_ms2(iT,1),phase(iT,1)+1);
         end
        
        % Make this calculation now in order to get Planet.Qmantle_Wm2 for making
        % filenames shortly       
        if CONVECTION_FLAG_I(iT) && eTBL_m(iT)>zb_outerIce_m(iT)
            warning('Convection predicted by not possible becuase the conductive layer thickness exceeds the thickness of the ice.')
            disp('Perhaps T_surf is outside the valid range for the scaling from Deschamps and Sotin 2001.')
            disp('Setting CONVECTION_FLAG_I to zero')
            CONVECTION_FLAG_I(iT) = 0;
        end
        % Convection has to be calculated prior to assigning depths in case
        % ice shell needs to be thinned to account for clathrates
   
        if CONVECTION_FLAG_I(iT)
            %conductive upper layer
            nConvectIce=(n_iceI(iT)+n_clath(iT))-nIceIIILithosphere-1; % indices where ice/claths exist in upper layer
            nconold=nConvectIce;
            if max_clath_depth<Planet.zb_outerIce_m(iT)
                %eTBL_m(iT)=eTBL_m(iT)+max_clath_depth; fixed in new
                %version
                T_K(iT,1:n_clath(iT))=linspace(Planet.Tsurf_K,T_clath_ice(iT),n_clath(iT));
                inds_eTBL = find(z_m(iT,n_clath(iT):nConvectIce)<=eTBL_m(iT));
                
                if ~isempty(inds_eTBL)
                    inds_eTBL=inds_eTBL +n_clath(iT);
                    if length(inds_eTBL)>1
                Pterm =(P_MPa(iT,inds_eTBL)-P_MPa(iT,inds_eTBL(1)))./(P_MPa(iT,inds_eTBL(end))-P_MPa(iT,inds_eTBL(1)));
           
                T_K(iT,inds_eTBL) = (Tc(iT).^(Pterm)).*(T_clath_ice(iT).^(1-Pterm));
                    else
                         T_K(iT,inds_eTBL)=T_clath_ice(iT)+Q_Wm2(iT)./kIce(iT).*(z_m(iT,inds_eTBL)-z_m(iT,inds_eTBL-1));
                    end
                else
                    inds_eTBL=n_clath(iT);
                end
            else
            inds_eTBL = find(z_m(iT,1:nConvectIce)<=eTBL_m(iT));
           
            Pterm = P_MPa(iT,inds_eTBL)./P_MPa(iT,inds_eTBL(end));
            if phase(iT,inds_eTBL(1))==1
            T_K(iT,inds_eTBL) = (Tc(iT).^(Pterm)).*(Planet.Tsurf_K.^(1-Pterm));
            elseif phase(iT,inds_eTBL(1))==30;
                 T_K(iT,1:inds_eTBL(end))=linspace(Planet.Tsurf_K,Tc(iT),length(inds_eTBL));
            end
            end
            
             rho_kgm3(iT,1:inds_eTBL(end))=getRhoIce(P_MPa(iT,1:inds_eTBL(end)),T_K(iT,1:inds_eTBL(end)),phase(iT,1:inds_eTBL(end)));
            
            P_bound=P_MPa(iT,inds_eTBL(end));
            T_K_ccbound=T_K(iT,inds_eTBL(end));
     
           
            %convective region

            for iconv = inds_eTBL(end)+1:nConvectIce % added parentheses around 1:nConvectIce 20200103 SDV

                rho_kgm3(iT,iconv)=getRhoIce(P_MPa(iT,iconv),T_K(iT,iconv-1),phase(iT,iconv));

                if POROUS_ICE % adjust if porosity needs to be considered
                    por_in(iT).p=P_MPa(iT,iconv)*1e-3;
                    por_in(iT).t = T_K(iT,iconv-1);
                    por_in(iT).den = rho;

                    if isfield(Planet,'phi_surface')
                        por_out(iT) = get_porosity_ice(por_in(iT),Planet.phi_surface);
                    else
                        por_out(iT) =get_porosity_ice(por_in(iT));
                    end
                    rho = por_out(iT).den;
                end
                    %[Cp(iT,iconv), alpha_K(iT,iconv)] = getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),1);
                    if phase(iT,iconv)==1 % if water ice, use SeaFreeze otherwise use Helgeurd.
                        cpout=SeaFreeze([P_MPa(iT,iconv),T_K(iT,iconv-1)],'Ih');
                        alpha_K(iT,iconv)=cpout.alpha;
                        Cp(iT,iconv)=cpout.Cp;
                    else

                        cpout=Helgerud_sI(P_MPa(iT,iconv),T_K(iT,iconv-1));
                        [Cp(iT,iconv) alpha_K(iT,iconv)]= getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),phase(iT,iconv)) ;
                        alpha_K(iT,iconv)=cpout.alpha; % use better alpha
                    end

                %aK = 1.56e-4; % thermal expansive? Switch and use SeaFreeze ( find clath values()

                T_K(iT,iconv) = T_K(iT,iconv-1)+alpha_K(iT,iconv).*T_K(iT,iconv)./Cp(iT,iconv)./rho_kgm3(iT,iconv)*deltaP*1e6;
                % double check temperatures make sense
                if strcmp(Planet.Ocean.comp,'NH3') % kluge svance july 7 2021. not sure why two phase tests are happening. generally, the seafreeze test should be performed when possible
                    phase_test = 1;
                else
                    phase_test = getIcePhase(P_MPa(iT,iconv),T_K(iT,iconv-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
                end
                phase_test2=SF_WhichPhase([P_MPa(iT,iconv),T_K(iT,iconv)]);

                if phase_test==0 || phase_test2==0 % & n_clath>0
                    Planet.Tb_K(iT)=T_K(iT,iconv);
                    %iconv=iconv-2;
                    Zm_break(iT)=z_m(iT,iconv); % saves the depth at which is broked
                    zb_outerIce_m(iT)=Zm_break(iT);   % new ice depth is saved

                    Zdiff(iT)=z_m(iT,nConvectIce)-Zm_break(iT); % change in depth
                     % previous convection layer

                    nConvectIce=iconv;% changes the indices for convection


                    %make adjustments to indices n_clath+n_ice+n_ocean
                    %should always be the same.
                    n_iceI(iT) = iconv - n_clath(iT); % new ice layer is the last index-n_clath
                    if n_iceI(iT)<0
                        n_iceI(iT) = 0; % indicates entire ice shell would be clathrates
                        iconv = iconv-1;
                        Zclath(iT) = z_m(iT,iconv);
                    end

                    n_ocean(iT) = nsteps - n_iceI(iT) - n_clath(iT);% adds more indices to ocean


                    PbI_MPa(iT)=P_MPa(iT,iconv);
                    Pb_MPa(iT)=P_MPa(iT,iconv);
                    %                        alculate parameters
%                     [Q_Wm2_new(iT),deltaTBL_m_new(iT),eTBL_m_new(iT),Tc_new(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),nu(iT),CONVECTION_FLAG_I(iT)]=...
%                         ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K(iT),PbI_MPa(iT)/2,zb_outerIce_m(iT),g_ms2(iT,1),phase(iT,iconv)+1);
                    inds_deltaTBL = find(z_m(iT,1:nConvectIce)>=z_m(iT,nConvectIce)-deltaTBL_m(iT));
                    T_K(iT,inds_deltaTBL) = (Planet.Tb_K(iT).^(P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT))).*(T_K(iT,inds_deltaTBL(1)-1).^(1-P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT)));

                    % %
                    %
                    deltaP = (Params.Pseafloor_MPa-Pb_MPa(iT))/n_ocean(iT); % in case PbMPa changed
                    % recalculate ocean
                    for il=1:n_ocean(iT)
                        ill=il+inds_deltaTBL(end);
                        % adjust for new ocean layers if necessary
                        P_MPa(iT,ill) = Pb_MPa(iT) + il*deltaP;
                        [rho_ocean,Cp(iT,ill),alpha_o]= fluidEOS(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);
                        if alpha_o<=0

                            disp('Ocean alpha at ice interface is less than zero. adjusting temperature upward.')
                            disp('This means there''s a conductive layer at the interface with thickness inversely proportional to the heat flow.')
                            disp('The thickness is likely less than a few 100 m. See Melosh et al. 2004.')
                            Tnew = fzero(@(T) alphaAdjust(P_MPa(iT,ill),T,wo,Planet.Ocean.comp),273);
                            disp(['The new temperature is ' num2str(Tnew) ', or delta T = ' num2str(Tnew-T_K(iT,ill-1)) '.'])
                            T_K(iT,ill)=Tnew;
                        else
                            T_K(iT,ill) = T_K(iT,ill-1)+ alpha_o*T_K(iT,ill-1)./(Cp(iT,ill))./(rho_ocean)*deltaP*1e6; %adiabatic gradient in ocean; this introduces an error, in that we are using the temperature from the previous step
                        end


                        phase(iT,ill) = getIcePhase(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI

                        if phase(iT,ill)==0
                            rho_kgm3(iT,ill) = rho_ocean;
                             [rho_ocean,Cp(iT,ill),alpha_o]= fluidEOS(P_MPa(iT,ill),T_K(iT,ill-1),wo,Planet.Ocean.comp);
                   
                        else
                            rho_kgm3(iT,ill) = getRhoIce(P_MPa(iT,ill),T_K(iT,ill),phase(iT,ill));
                            
                            if ~isfield(Planet.Ocean,'fnTfreeze_K')

                                T_K(iT,ill) = getTfreeze(P_MPa(iT,ill),wo,Planet.Ocean.comp,T_K(iT,ill-1)); %find the temperature on the liquidus line corresponding to P(k,i+nIceI); should really use a conductive profile here, but that would seem to naturally bring us back to the liquidus. Steve wants to confirm.
                            else
                                T_K(iT,ill) = Planet.Ocean.fnTfreeze_K(P_MPa(iT,ill),wo);
                            end
                            if T_K(iT,ill)<T_K(iT,ill-1)
                                T_K(iT,ill) = T_K(iT,ill-1);
                            end
                           [Cp(iT,ill) alpha_K(iT,ill)]= getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),phase(iT,ill)) ;
                  

                        end

                    end
                    z_ocean_m=z_m(iT,ill);

                    break
                end
            end
            
            if nconold==nConvectIce
                % bottom layer of conduction, recalculates for phase at bottom
                 [Q_Wm2(iT),deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),nu(iT),CONVECTION_FLAG_I(iT)]=...
                        ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K(iT),PbI_MPa(iT)/2,zb_outerIce_m(iT),g_ms2(iT,1),phase(iT,iconv)+1);

                inds_deltaTBL = find(z_m(iT,1:nConvectIce)>=z_m(iT,nConvectIce)-deltaTBL_m(iT));
                %
                %
                T_K(iT,inds_deltaTBL) = (Planet.Tb_K(iT).^(P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT))).*(T_K(iT,inds_deltaTBL(1)-1).^(1-P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT)));
                %           end
                rho_kgm3(iT,inds_deltaTBL) = getRhoIce(P_MPa(iT,inds_deltaTBL),T_K(iT,inds_deltaTBL),phase(iT,inds_deltaTBL));
                 %[Cp(iT,inds_deltaTBL) alpha_K(iT,inds_deltaTBL)]= getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),phase(iT,inds_deltaTBL)) ;
                        

                z_ocean_m(iT)= z_m(iT,inds_deltaTBL(end)+1);

                if POROUS_ICE
                    por_in(iT).p=P_MPa(iT,inds_deltaTBL)*1e-3;
                    por_in(iT).t = T_K(iT,inds_deltaTBL);
                    por_in(iT).den = rho_kgm3(iT,inds_deltaTBL);
                    if isfield(Planet,'phi_surface')
                        por_out(iT) = get_porosity_ice(por_in(iT),Planet.phi_surface);
                    else
                        por_out(iT) =get_porosity_ice(por_in(iT));
                    end

                    rho_kgm3(iT,inds_deltaTBL) = por_out(iT).den;
                end
            end
        else
            z_ocean_m(iT) = zb_outerIce_m(iT);

%             if nconold==nConvectIce
%                 % bottom layer of conduction, recalculates for phase at bottom
%                  [Q_Wm2(iT),deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),nu(iT),CONVECTION_FLAG_I(iT)]=...
%                         ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K(iT),PbI_MPa(iT)/2,zb_outerIce_m(iT),g_ms2(iT,1),phase(iT,iconv)+1);
% 
%                 inds_deltaTBL = find(z_m(iT,1:nConvectIce)>=z_m(iT,nConvectIce)-deltaTBL_m(iT));
%                 %
%                 %
%                 T_K(iT,inds_deltaTBL) = (Planet.Tb_K(iT).^(P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT))).*(T_K(iT,inds_deltaTBL(1)-1).^(1-P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT)));
%                 %           end
%                 rho_kgm3(iT,inds_deltaTBL) = getRhoIce(P_MPa(iT,inds_deltaTBL),T_K(iT,inds_deltaTBL),phase(iT,inds_deltaTBL));
%                  %[Cp(iT,inds_deltaTBL) alpha_K(iT,inds_deltaTBL)]= getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),phase(iT,inds_deltaTBL)) ;
%                         
% 
%                 z_ocean_m(iT)= z_m(iT,inds_deltaTBL(end)+1);
% 
%                 if POROUS_ICE
%                     por_in(iT).p=P_MPa(iT,inds_deltaTBL)*1e-3;
%                     por_in(iT).t = T_K(iT,inds_deltaTBL);
%                     por_in(iT).den = rho_kgm3(iT,inds_deltaTBL);
%                     if isfield(Planet,'phi_surface')
%                         por_out(iT) = get_porosity_ice(por_in(iT),Planet.phi_surface);
%                     else
%                         por_out(iT) =get_porosity_ice(por_in(iT));
%                     end
% 
%                     rho_kgm3(iT,inds_deltaTBL) = por_out(iT).den;
%                 end
%             end

        end % if CONVECTION_FLAG
        
        
        if Planet.EQUIL_Q
            if (~isfield(Params,'NO_ICEI_CONVECTION') || Params.NO_ICEI_CONVECTION == false) && ...
               (CONVECTION_FLAG_I(iT) || (isfield(Params,'FORCE_ICEI_CONVECTION') && Params.FORCE_ICEI_CONVECTION == true))
                Qmantle_Wm2(iT) = Q_Wm2(iT);
            else
                Qmantle_Wm2(iT) = Qb(iT);
            end
        else
            Qmantle_Wm2(iT) = 2.2e11/4/pi/Planet.R_m^2; % this is more reasonable for radiogenic only
        end
        Planet.Qmantle_Wm2(iT) = Qmantle_Wm2(iT);
        
        deltaP = (Params.Pseafloor_MPa-Pb_MPa(iT))/n_ocean(iT); %
        
        Planet.zClath_m(iT) = Zclath(iT);
     
        Planet.Profile_ID(iT) = ['Ts' num2str(Planet.Tsurf_K,'%0.0f') 'Zb' strLow num2str(Planet.zb_outerIce_m(iT),'%0.0f') ...
            'mQm' num2str(1000*Planet.Qmantle_Wm2(iT),'%0.0f') 'mWm2_CMR2p' ...
            num2str(10000*Planet.Cmeasured,'%0.0f') '_' thiseos];
        Planet.Profile_fname(iT) = [savefile '_' char(Planet.Profile_ID(iT))];
            
            savestr = [savebase Planet.Ocean.comp ...
                '_' num2str(10*round(Planet.Ocean.w_ocean_pct)) 'WtPpt' minEOS porIceStr ...
                '_Tb' num2str(round(Planet.Tb_K(iT),3),'%.3f') 'K' ];
            saveFile = fullfile([datpath savestr '.txt']);
            %save(fullfile([datpath savefile '_pp' num2str(iT)]),'P_MPa','Pb_MPa','PbI_MPa','nIceIIILitho','T_K','Tb_K','phase','deltaP','wo','nTbs','rho_kgm3','rho_ocean','Cp','alpha_o','nsteps','n_clath','n_iceI','n_ocean','max_clath_depth'); % save the progress at each step
            dlmwrite(saveFile, '  nHeadLines = 15', 'delimiter', '');
            dlmwrite(saveFile, ['  Tb_K = ' num2str(Planet.Tb_K(iT))], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  Zb_km = ' num2str(Planet.zb_outerIce_m(iT)/1e3)], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  zClath_m = ' num2str(Planet.zClath_m(iT))], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  Pb_MPa = ' num2str(Pb_MPa(iT))], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  PbI_MPa = ' num2str(PbI_MPa(iT))], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  deltaP = ' num2str(deltaP)], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  alpha_o = ' num2str(alpha_o)], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  Qmantle_Wm2 = ' num2str(Planet.Qmantle_Wm2(iT))], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  nStepsIceI = ' num2str(Params.nsteps_iceI)], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  nStepsOcean = ' num2str(Params.nsteps_ocean)], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  nStepsHydro = ' num2str(nsteps)], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  nIceIIILitho = ' num2str(nIceIIILithosphere)], 'delimiter', '', '-append');
            dlmwrite(saveFile, ['  nIceVLitho = ' num2str(nIceVLithosphere)], 'delimiter', '', '-append');
            header = sprintf('%s\t\t%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\t\t',...
            'z (km)', 'r (m)', 'P (MPa)', 'T (K)', 'phase ID', 'rho (kg/m3)', 'g (m/s2)', 'Cp (W/kg/K)');
            dlmwrite(saveFile, header, 'delimiter', '', '-append');

            saveData = [z_m(iT,:)'*1e-3 r_m(iT,:)' P_MPa(iT,:)' T_K(iT,:)' phase(iT,:)' rho_kgm3(iT,:)' g_ms2(iT,:)' Cp(iT,:)'];
            dlmwrite(saveFile,saveData,...
                'delimiter','\t',...
                'precision',16,...
                '-append');
        end
    else
        [z_m, r_m, P_MPa, T_K, phase, rho_kgm3, g_ms2, Cp, rho_ocean] = deal(zeros(nTbs,nsteps));
        for iT = 1:nTbs
            
            savestr = [savebase Planet.Ocean.comp ...
                '_' num2str(10*round(Planet.Ocean.w_ocean_pct)) 'WtPpt' minEOS porIceStr ...
                '_Tb' num2str(round(Planet.Tb_K(iT),3),'%.3f') 'K' ];
            saveFile = fullfile([datpath savestr '.txt']);
            %try
                %load(fullfile([datpath savefile '_pp' num2str(iT)]));
                
                % Parse header from text file
                fReload = fopen(saveFile);
                    nHeadLinesStr = split(fgetl(fReload),'=');
                    Tb_KStr = split(fgetl(fReload),'=');
                    Zb_kmStr = split(fgetl(fReload),'=');
                    zClath_mStr = split(fgetl(fReload),'=');
                    Pb_MPaStr = split(fgetl(fReload),'=');
                    PbI_MPaStr = split(fgetl(fReload),'=');
                    deltaPStr = split(fgetl(fReload),'=');
                    alpha_oStr = split(fgetl(fReload),'=');
                    Qmantle_Wm2Str = split(fgetl(fReload),'=');
                    nStepsIceIStr = split(fgetl(fReload),'=');
                    nStepsOceanStr = split(fgetl(fReload),'=');
                    nStepsHydroStr = split(fgetl(fReload),'=');
                    nIceIIILithoStr = split(fgetl(fReload),'=');
                    nIceVLithoStr = split(fgetl(fReload),'=');
                    
                    nHeadLines = str2num(nHeadLinesStr{2});
                    Planet.Tb_K(iT) = str2num(Tb_KStr{2});
                    Planet.zb_outerIce_m(iT) = str2num(Zb_kmStr{2});
                    Planet.zClath_m(iT) = str2num(zClath_mStr{2});
                    Pb_MPa = str2num(Pb_MPaStr{2});
                    PbI_MPa = str2num(PbI_MPaStr{2});
                    deltaP = str2num(deltaPStr{2});
                    alpha_o = str2num(alpha_oStr{2});
                    Planet.Qmantle_Wm2(iT) = str2num(Qmantle_Wm2Str{2});
                    nStepsIceI = str2num(nStepsIceIStr{2});
                    nStepsOcean = str2num(nStepsOceanStr{2});
                    nStepsHydro = str2num(nStepsHydroStr{2});
                    nIceIIILitho = str2num(nIceIIILithoStr{2});
                    nIceVLitho = str2num(nIceVLithoStr{2});
                fclose(fReload);
                
                reloadData = dlmread(saveFile, '', nHeadLines);
                if ~exist('max_clath_depth')
                    max_clath_depth = 1e15;
                end
            %catch
                %error(['ERROR: cfg.CALC_NEW=0 but ' savefile ' was not found. Re-run with cfg.CALC_NEW set to 1 to generate the Profile.']);
            %end
           
            z_m(iT,:) = 1e3 * reloadData(:,1);
            r_m(iT,:) = reloadData(:,2);
            P_MPa(iT,:) = reloadData(:,3);
            T_K(iT,:) = reloadData(:,4);
            phase(iT,:) = reloadData(:,5);
            rho_kgm3(iT,:) = reloadData(:,6);
            g_ms2(iT,:) = reloadData(:,7);
            Cp(iT,:) = reloadData(:,8);
            
            %'nsteps','n_clath','n_iceI','n_ocean'
            nIceIIILithosphere = nIceIIILitho;
            nIceVLithosphere = nIceVLitho;
            n_clath = nStepsHydro - nStepsIceI - nStepsOcean;
            n_iceI(iT) = nStepsIceI;
            n_ocean = nStepsOcean;
            rho_ocean(iT,nStepsIceI+1:nStepsIceI+nStepsOcean) = rho_kgm3(iT,nStepsIceI+1:nStepsIceI+nStepsOcean);
            Tb_K = Planet.Tb_K(iT);
            save(fullfile([datpath savefile '_pp' num2str(iT)]),'P_MPa','Pb_MPa','PbI_MPa','nIceIIILithosphere','T_K','Tb_K','phase','deltaP','wo','nTbs','rho_kgm3','rho_ocean','Cp','alpha_o','nsteps','n_clath','n_iceI','n_ocean','max_clath_depth');
        end
    end
    
end




% Assign electrical conductivities for the ocean
if cfg.CONDUCT
    k_S_m = NaN*ones(size(P_MPa));
    switch Planet.Ocean.comp
    case 'MgSO4' 
        if wo==10 && isfield(Planet.Ocean,'Electrical') && strcmp(Planet.Ocean.Electrical,'Pan2020')
        % use values from Pan et al. 2020
        for iTb = 1:length(Planet.Tb_K)
            thisP = P_MPa(iTb,phase(iTb,:)==0);
            thisT = T_K(iTb,phase(iTb,:)==0);
            for iP = 1:length(thisP)
               kvals(iP)  = getSigmaMgSO4_Pan(thisP(iP),thisT(iP));
            end
                k_S_m(iTb,phase(iTb,:)==0) = kvals;
        end
        else
        % Extrapolation from Larionov 1984 to high concentration and low
        % temperature
        LarionovMgSO4 = getLarionov1984MgSO4Conductivities(1);
        Pextrap = LarionovMgSO4.Pextrap_MPa;
        Textrap = LarionovMgSO4.Textrap_K;
        k_S_m_extrap = LarionovMgSO4.k_S_m_extrap_p01m;
        [Pgg,Tgg]=meshgrid(Pextrap,Textrap);
        for iTb = 1:length(Planet.Tb_K)
            k_S_m(iTb,phase(iTb,:)==0) = interp2(Pextrap,Textrap,k_S_m_extrap,P_MPa(iTb,phase(iTb,:)==0),T_K(iTb,phase(iTb,:)==0),'spline') * conduct_scaling_MgSO4;
        end
        end
    case 'Seawater'
        if wo>2
            for iTb = 1:length(Planet.Tb_K)
                k_S_m(iTb,phase(iTb,:)==0)=swEOS.gsw.C(wo*ones(1,length(T_K(iTb,phase(iTb,:)==0))),T_K(iTb,phase(iTb,:)==0),P_MPa(iTb,phase(iTb,:)==0)*10);
            end
        end
    end
end


%% Deeper Interior
MR2 = Planet.M_kg*Planet.R_m^2;
dz_m = gradient(z_m);
nz = length(dz_m);
nR = nz;
% nR = round(3/4*nz);
npre = nz-nR; % confine the search for structure to the lower part
 C_H2O = sum((8/3*pi*rho_kgm3(:,1:npre).*r_m(:,1:npre).^4.*dz_m(:,1:npre))')'; %calculate C_H2O to the beginning of the depth where we start probing for the beginning of the silicate layer

% Preallocate
C2inds = cell(1,nTbs);
[C2mean, C2max, C2min, R_sil_mean_m, M_H2O_mean_kg] = deal(zeros(1,nTbs));
[R_sil_m, rho_sil_kgm3] = deal(zeros(nTbs,nR));
    
    % without a core
if ~Planet.FeCore
    C1 = zeros(nTbs,nR);
    for iz = 1:nR
        C_H2O(:,iz+1) = C_H2O(:,iz)+(8/3*pi*rho_kgm3(:,iz+npre).*r_m(:,iz+npre).^4.*dz_m(:,iz+npre));
        R_sil_m(:,iz) = r_m(:,iz+npre);
        rho_sil_kgm3(:,iz) = 3/4/pi*(Planet.M_kg-M_above_kg(:,iz+npre))./power(R_sil_m(:,iz),3);
        C1(:,iz) = C_H2O(:,iz+1)+8/15*pi*power(R_sil_m(:,iz),5).*rho_sil_kgm3(:,iz);
    end
    for iT = 1:nTbs
        C2inds{iT} = find(C1(iT,:)/MR2>Planet.Cmeasured-Planet.Cuncertainty & C1(iT,:)/MR2<Planet.Cmeasured+Planet.Cuncertainty);
        if isempty(C2inds{iT})
            error(['C/MR2=' num2str(Planet.Cmeasured) ' not found. min:' num2str(min(C1(iT,:)/MR2)) '; max:' num2str(max(C1(iT,:)/MR2))])
        end
        C2mean(iT) = round(mean(C2inds{iT}));
        C2max(iT) = max(C2inds{iT});
        C2min(iT) = min(C2inds{iT});
        R_sil_mean_m(iT) = R_sil_m(iT,C2mean(iT));
        M_H2O_mean_kg(iT) = M_above_kg(iT,C2mean(iT));
    end
    R_Fe_mean_m = zeros(1,nTbs);
   
%     dz_ocean_m_m(~dindsVI & ~dindsV) = Planet.R_m - R_sil_mean_m(~dindsVI &
%     ~dindsV)-zI_m(~dindsVI & ~dindsV); not sure this is right. commenting
%     out on March 8 2018

    if ~cfg.SKIP_PROFILES
        % Plot the results
        set(0, 'CurrentFigure', figs.mant);
        clf;hold all
        set(gcf,'Position', [335 133 696 547], 'Name', lbl.mantl)
        for iT=1:nTbs
            plot(rho_sil_kgm3(iT,C2inds{iT})',R_sil_m(iT,C2inds{iT})'*1e-3);
        end
        lstr_3 = cell(1,nTbs);
        for iT = 1:nTbs
            lstr_3{iT} = [math 'T_{b}' nm ': ' num2str(Planet.Tb_K(iT),'%0.1f') ' K'];
        end
        legend(lstr_3,'Fontsize',lbl.smtext)
        box on
        xlabel([math '\rho_{' nm 'sil}' nm ' (kg m^{-3})'],'Fontsize',lbl.mltext);
        ylabel([math 'R_{' nm 'sil}' nm ' (km)'],'Fontsize',lbl.mltext)
        title(['No Fe core ; ' math 'C/MR^2' bnm ' = ' num2str(Planet.Cmeasured) '\pm' num2str(Planet.Cuncertainty) ';' math ' W ' bnm ' = ' num2str(wo) ' wt%'],'Fontsize',lbl.lgtext)

        print(figs.mant,Params.savefigformat,fullfile([figpath savebase vmant cfg.xtn]));
    end
% =====
else % WITH A CORE
    rho_Fe_kgm3 = Planet.rhoFe*Planet.rhoFeS/(Planet.xFeS*(Planet.rhoFe-Planet.rhoFeS)+Planet.rhoFeS);
    [C2,M_iron_kg,R_Fe_m] = deal(zeros(nTbs,nR));
    for iz = 1:nR
        C_H2O(:,iz+1) = C_H2O(:,iz)+(8/3*pi*rho_kgm3(:,iz+npre).*r_m(:,iz+npre).^4.*dz_m(:,iz+npre));
        R_sil_m(:,iz) = r_m(:,iz+npre);
        for iT = 1:nTbs
            [C2(iT,iz),R_Fe_m(iT,iz)] = CoreSize(Planet.rho_sil_withcore_kgm3,rho_Fe_kgm3,C_H2O(iT,iz+1),M_above_kg(iT,iz+npre),R_sil_m(iT,iz), Planet.M_kg, Planet.R_m);
        end
    end
    
    [R_Fe_mean_m, R_Fe_range_m, R_sil_range_m] = deal(zeros(1,nTbs));
    r_core_m = zeros(nTbs,Params.nsteps_core);
    lstr_3 = cell(1,nTbs);
    for iT=1:nTbs
        C2inds{iT} = find(C2(iT,:)/MR2>Planet.Cmeasured-Planet.Cuncertainty & C2(iT,:)/MR2<Planet.Cmeasured+Planet.Cuncertainty);
        if isempty(C2inds{iT})
            error(['C/MR2=' num2str(Planet.Cmeasured) ' not found. min:' num2str(min(C2(iT,:)/MR2)) '; max:' num2str(max(C2(iT,:)/MR2))])
        end
        C2mean(iT) = round(mean(C2inds{iT}));
        C2max(iT) = max(C2inds{iT});
        C2min(iT) = min(C2inds{iT});
        R_Fe_mean_m(iT) = R_Fe_m(iT,C2mean(iT));
        R_sil_mean_m(iT) = R_sil_m(iT,C2mean(iT));
        M_H2O_mean_kg(iT) = M_above_kg(iT,C2mean(iT));
        R_Fe_range_m(iT) = R_Fe_m(iT,C2max(iT))-R_Fe_m(iT,C2min(iT));
        R_sil_range_m(iT) = R_sil_m(iT,C2min(iT))-R_sil_m(iT,C2max(iT));
        
        r_core_m(iT,:) = linspace(R_Fe_mean_m(iT),0,Params.nsteps_core);
        if ~cfg.SKIP_PROFILES; lstr_3{iT} = [math 'T_{b}' nm ': ' num2str(Planet.Tb_K(iT),'%0.1f') ' K']; end
    end

    if ~cfg.SKIP_PROFILES
        set(0, 'CurrentFigure', figs.core);
        clf;hold all
        set(gcf,'Position', [335 133 854 547], 'Name', lbl.corsz)
        
        plot(R_Fe_m(iT,C2inds{iT})'*1e-3,R_sil_m(iT,C2inds{iT})'*1e-3);
        legend(lstr_3,'Fontsize',lbl.smtext)
        box on
        xlabel([math 'R_{' nm 'Fe}' nm ' (km)'],'Fontsize',lbl.mltext);
        ylabel([math 'R_{' nm 'sil}' nm ' (km)'],'Fontsize',lbl.mltext);
        title(['Fe core ; ' math 'C/MR^2' bnm ' = ' num2str(Planet.Cmeasured) '\pm' num2str(Planet.Cuncertainty) '; ' math ' W ' bnm ' = ' num2str(wo) ' wt%; ' math '\rho_{' bnm 'sil}' bnm ': ' num2str(Planet.rho_sil_withcore_kgm3,'%0.0f') '; ' math '\rho_{' bnm 'Fe}' bnm ': ' num2str(rho_Fe_kgm3,'%0.0f')],'Fontsize',lbl.lgtext)

        print(figs.core,Params.savefigformat,fullfile([figpath savebase vcore cfg.xtn]));
    end
end


%% Print the depths to the various layers
%allocate
zI_m = Planet.zb_outerIce_m;
%zI_m = zb_outerIce_m-Zclath;
[zIII_m, zV_m, zVI_m, indSil] = deal(zeros(1,nTbs));
indsV = zV_m;
indsVI = zV_m;
indsIII = zV_m;
indsClath = zV_m;
dz_ocean_m_m = Planet.zb_outerIce_m;

    % figure out the indices for the tops of the different ice layers
for iT = 1:nTbs
      theseinds = find(phase(iT,:)==1);
    if ~isempty(theseinds)
        indsI(iT) = theseinds(1);
        zI_m2(iT) = z_m(iT,indsI(iT));
    elseif isempty(theseinds) & Params.nsteps_clath>0
        %indsI(iT) = theseinds(1);
        zI_m2(iT) = 0;
    end
    
    theseinds = find(phase(iT,:)==3);
    if ~isempty(theseinds)
        indsIII(iT) = theseinds(1);
        zIII_m(iT) = z_m(iT,indsIII(iT));
    end
    theseinds = find(phase(iT,:)==5);
    if ~isempty(theseinds)
        indsV(iT) = theseinds(1);
        zV_m(iT) = z_m(iT,indsV(iT));
    end
    theseinds = find(phase(iT,:)==6);
    if ~isempty(theseinds)
        indsVI(iT) = theseinds(1);
        zVI_m(iT) = z_m(iT,indsVI(iT));
    end
    
    indSil(iT) = find(R_sil_m(iT,:)==R_sil_mean_m(iT));
end

% This step makes sure that the mantle length matrices are all the same
% length, by setting indSil(iT) + nsteps_mantle(iT) = Params.nsteps_mantle.
if nTbs == 1
    nsteps_mantle = Params.nsteps_mantle;
else
    nsteps_mantle = [Params.nsteps_mantle Params.nsteps_mantle*ones(1,nTbs-1)+(indSil(1)-indSil(2:end))];%
    if find(nsteps_mantle<0); error('Problems indexing in the "mantle" layer. Consider increasing nsteps_mantle.'); end
end

[dz_ocean_m_m,dzIII_m,dzV_m,dzVI_m,dzClath_m] = deal(zeros(1,nTbs));
%find the radii at the tops of the different layers
RIII_m =Planet.R_m-zIII_m;
RV_m = Planet.R_m-zV_m;
RVI_m = Planet.R_m-zVI_m;

Rclath_m=Planet.R_m-Planet.zClath_m;
Rice_m=Planet.R_m-zI_m2;

% find the thicknesses of the layers.  keep adjusting the thicknesses
% accordingly
dindsclath = Rclath_m>R_sil_mean_m & Planet.zClath_m>0;
dzClath_m(dindsclath) = Planet.zClath_m(dindsclath)-zI_m(dindsclath);


dindsIII = RIII_m>R_sil_mean_m & zIII_m>0;
dz_ocean_m_m(dindsIII) = zIII_m(dindsIII)-zI_m(dindsIII);

dindsV = RV_m>R_sil_mean_m & zV_m>0;
dz_ocean_m_m(~dindsIII & dindsV) = zV_m(~dindsIII & dindsV)-zI_m(~dindsIII & dindsV);
dzIII_m(dindsIII & dindsV) = zV_m(dindsV & dindsIII)- zIII_m(dindsV & dindsIII);

dindsVI = RVI_m>R_sil_mean_m & zVI_m>0;
dz_ocean_m_m(dindsVI & ~dindsV) = zVI_m(dindsVI & ~dindsV) - zI_m(dindsVI & ~dindsV);

% use R_sil_mean_m to calculate thickness of ice VI
for iT = 1:nTbs
    if ~isempty(dindsVI(iT)) && dindsVI(iT)>0
        dzVI_m(iT) = Planet.R_m-R_sil_mean_m(iT) - zVI_m(iT);
    elseif ~isempty(dindsV(iT)) && dindsV(iT)>0
        dzVI_m(iT) = 0;
        dzV_m(iT) = RV_m(iT)-RVI_m(iT)*dindsVI(iT); % this "fix" may introduce an failure condition, but the previous method on the next line was also failing.
%         dzV_m(iT) = RV_m(iT)-R_sil_mean_m(iT) - zVI_m(iT)*dindsVI(iT);
    elseif ~isempty(dindsIII(iT)) && dindsIII(iT)>0
        dzV_m(iT) = 0;
        dzVI_m(iT) = 0;
        dzIII_m(iT) = Planet.R_m-R_sil_mean_m(iT) - zVI_m(iT);
%   elseif ~isempty(dindsII(iT)) && dindsII(iT)>0
%       dzV_m(iT) = 0;
%       dzVI_m(iT) = 0;
%       dzIII_m(iT) = Planet.R_m-R_sil_mean_m(iT) - zVI_m(iT);
    end
end

dz_ocean_m_m(~dindsVI & ~dindsV) = Planet.R_m - R_sil_mean_m(~dindsVI & ~dindsV)-zI_m(~dindsVI & ~dindsV);
zTotal_m = zI_m+dz_ocean_m_m+dzIII_m+dzV_m+dzVI_m+abs(dzClath_m);


    % this should be elaborated upon to compute the actual mass by
    % integrating the masses of the concentric shells and removing the
    % effect of the salt. This can be easily done using eqst. It will be
    % interesting to compare the mass of the pure water ice and ocean layer
    % with that of the salty layer
    % less precise:
%     MH2O = 4/3*pi*(Planet.R_m.^3 - R_sil_mean_m.^3).*mean(rho_kgm3(:,C2mean));
    % more precise:
WtH2O = M_H2O_mean_kg/Planet.M_kg; % this is the predicted water in the planet
dWtH2O = Planet.XH2O-WtH2O;
    
if Planet.FeCore
    Mcore = 4/3*pi*R_Fe_mean_m.^3.*rho_Fe_kgm3;
    Wtcore = Mcore/Planet.M_kg;
    dWtcore = Planet.xFe_core-Wtcore;
    XS = 32.065/87.91*Planet.xFeS; %mass fraction of S in the core
    Msulfur = Mcore.*XS; % mass of sulfur in the core
    WtS = Msulfur/Planet.M_kg; % as a fraction of Europa's mass
    dWtS = Planet.xFeS_meteoritic*32.065/87.91-WtS;
end

% We now have everything we need for LayeredInductionResponse, iff we are
% setting interior conductivities naively. If we care about modeling detailed
% mantle and core conductivities, we need to get further along in calculating the depth
% profile (including mantle and core properties) in order to assign those values.
nsteps_tot = Params.nsteps_core + nsteps_mantle(1) + indSil(1)-1;
iceSig = 1e-16;
mantleSig = 1e-16;
if ~isfield(Planet, 'coreSig')
    coreSig = 1e-16;
else
    coreSig = Planet.coreSig;
end

if ~Planet.FeCore; nsteps_tot = nsteps_tot - Params.nsteps_core; end
[r_Planet_m, sig] = deal(zeros(nTbs,nsteps_tot));
[ocean_thk, ind_Ih, ind_Obot, kmean, ktop] = deal(zeros(1,nTbs));
Planet.ice_thk = strings(1,nTbs);
totL = size(r_Planet_m,2);
for iT=1:nTbs
    % Grab the indices for the ices, but only for the portion of the
    % interior calculation above the silicate interface
    H2Oinds = 1:indSil(iT)-1;
    indsI = find(phase(iT,H2Oinds)==1);
    indsSurf = find(phase(iT,H2Oinds)==1 | phase(iT,H2Oinds)==30);% changes to look for top ice/clath layer
    indsClath = find(phase(iT,H2Oinds)==30);
    indsLiquid = find(phase(iT,H2Oinds)==0);
    indsIII = find(phase(iT,H2Oinds)==3);
    indsV = find(phase(iT,H2Oinds)==5);
    indsVI = find(phase(iT,H2Oinds)==6);
    sig(iT,indsI) = iceSig;
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        sig(iT,indsLiquid) = k_S_m(iT,indsLiquid);
    else
        sig(iT,indsLiquid) = 0.0;
    end
    sig(iT,indsIII) = iceSig;
    sig(iT,indsV) = iceSig;
    sig(iT,indsVI) = iceSig;
    sig(iT,indSil(iT):indSil(iT)+nsteps_mantle(iT)-1) = mantleSig;
    if Planet.FeCore
        R_sil_bot = R_Fe_mean_m(iT) - (R_Fe_mean_m(iT) - R_sil_mean_m(iT))/(nsteps_mantle(iT)-1);
        r_mantle_m = linspace(R_sil_mean_m(iT),R_sil_bot,nsteps_mantle(iT));
        sig(iT,(nsteps_tot - Params.nsteps_core):end) = coreSig;
        r_Planet_m(iT,:) = [r_m(iT,1:indSil(iT)-1) r_mantle_m r_core_m(iT,:)];
    else
        r_mantle_m = linspace(R_sil_mean_m(iT),0,nsteps_mantle(iT));
        r_Planet_m(iT,:) = [r_m(iT,1:indSil(iT)-1) r_mantle_m];
    end
    sig(sig==0) = 1e-16; % zeros make the integration unstable
    sig(isnan(sig)) = 1e-16;

    ocean_thk(iT) = r_Planet_m(iT,indsLiquid(1)) - r_Planet_m(iT,indsLiquid(end)+1);
    disp(['Ocean thickness: ' num2str(ocean_thk(iT)*1e-3,'%0.1f km')])
    if ~isempty(indsI)
        ind_Ih(iT) = indsI(end);
    else
        ind_Ih(iT) = indsSurf(end);
    end
    ind_Obot(iT) = indsLiquid(end);
    kmean(iT) = mean(sig(iT,ind_Ih(iT)+1:ind_Obot(iT)));
    ktop(iT) = sig(ind_Ih(iT)+1);
    
    
    %% Interpolate fewer ocean layers to reduce computational load
    if cfg.REDUCED && Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        nIntL = cfg.nIntL;
        ocStart = totL - indsLiquid(end) + 1;
        ocEnd = totL - indsLiquid(1) + 1;
        nOis = length(indsLiquid);
        interpBds = zeros(1,length(r_Planet_m(iT,:)) - nOis + nIntL);
        interpSig = zeros(1,length(r_Planet_m(iT,:)) - nOis + nIntL);

        r_asc = flip(r_Planet_m(iT,:),2);
        sig_asc = flip(sig(iT,:),2);

        % Preserve non-ocean conductivities as-is
        interpBds(1:ocStart) = r_asc(1:ocStart);
        interpBds(ocStart+nIntL:end) = r_asc(ocEnd+1:end);
        interpSig(1:ocStart-1) = sig_asc(1:ocStart-1);
        interpSig(ocStart+nIntL:end) = sig_asc(ocEnd+1:end);
        if nOis == ocEnd - ocStart + 1
            % Single ocean layer, simple layer averages
            fullOcean = linspace(r_asc(ocStart),r_asc(ocEnd+1),nIntL+1);
            interpBds(ocStart:ocStart+nIntL-1) = fullOcean(2:end);
            fullSigs = zeros(1,nIntL);
            for iR = 1:nIntL
                avg_start = floor((iR-1)*nOis/nIntL);
                avg_end = floor(iR*nOis/nIntL) - 1;
                fullSigs(iR) = mean(sig_asc(ocStart+avg_start:ocStart+avg_end));
            end
            interpSig(ocStart:ocStart+nIntL-1) = fullSigs;
        else
            % Multiple ocean layers
            warning('WARNING: Multiple ocean layers found. Interpolate by hand.')
        end

        % Now pad interpolated arrays to make Planet fields consistent lengths
        padL = totL - length(interpBds);
        padSig = zeros(1,padL) + 1e-16;
        padBds = linspace(interpBds(end)-(padL+1)*1e-5, interpBds(end)-1e-5, padL);
        Planet.boundaries(iT,:) = [interpBds(1:end-1) padBds interpBds(end)];
        Planet.sig(iT,:) = [interpSig(1:end-1) padSig interpSig(end)];
    else
        Planet.boundaries(iT,:) = flip(r_Planet_m(iT,:),2);
        Planet.sig(iT,:) = flip(sig(iT,:),2);
    end
end

Planet.kmean(iT,:) = kmean(:);
Planet.ktop(iT,:) = ktop(:);

Planet.D_Ih_km = Planet.zb_outerIce_m(:)/1e3;
Planet.D_ocean_km = ocean_thk(:)/1e3;
Planet.ice_thk = string(num2str(Planet.D_Ih_km(:),'%0.0f'));
Planet.salt = string([Planet.Ocean.comp ' ' num2str(Planet.Ocean.w_ocean_pct) ' wt%']);
Planet.Ocean.indTop = nsteps_tot - ind_Ih(:); % Index in Planet.boundaries for outer radius of ocean (ice-ocean boundary)
Planet.Ocean.indSil = nsteps_tot - ind_Obot(:); % Index in Planet.boundaries for outer radius of silicate mantle
Planet.Ocean.indBot = nsteps_tot - ind_Obot(:) + 1; % Index in Planet.boundaries for outer radius of first ocean layer

% Layered induction not yet implemented for other than the Galilean moons.
if (strcmp(Planet.name,'Europa') || strcmp(Planet.name,'Ganymede') || strcmp(Planet.name,'Callisto'))  && isfield(Planet,'peaks_Hz')
    Planet.peaks_hr = 1./Planet.peaks_Hz/3600;
end
    
% Exit now if we intend to skip creating profiles and just do layered
% induction calculations.
if cfg.SKIP_PROFILES && ~cfg.CALC_NEW
    outPlanet = Planet;
    return
end

% Print layer depths for user
if cfg.DISP_LAYERS
    disp(['Tb:                    ' num2str(Planet.Tb_K,'\t%0.2f')])
    disp(['z(km) ice I:           ' num2str(zI_m*1e-3,'\t%0.1f')])
    disp(['z(km) clath:           ' num2str(Planet.zClath_m*1e-3,'\t%0.1f')])
    disp(['z(km) ice III:         ' num2str(dindsIII.*zIII_m*1e-3,'\t%0.0f')])
    disp(['z(km) ice V:           ' num2str(dindsV.*zV_m*1e-3,'\t%0.0f')])
    disp(['z(km) ice VI:          ' num2str(dindsVI.*zVI_m*1e-3,'\t%0.0f')])
    disp(['dz(km) Ocean:          ' num2str(dz_ocean_m_m*1e-3,'\t%0.0f')])
    disp(['dz(km) ice III:        ' num2str(dzIII_m*1e-3,'\t%0.0f')])
    disp(['dz(km) ice V:          ' num2str(dzV_m*1e-3,'\t%0.0f')])
    disp(['dz(km) ice VI:         ' num2str(dzVI_m*1e-3,'\t%0.0f')])
    disp(['dz(km) ice V + ice VI: ' num2str((dzV_m+dzVI_m)*1e-3,'\t%0.0f')])
    if isfield(Planet,'XH2O')
    disp('')
    disp(['wt Pct of H2O:         ' num2str(100*WtH2O,'\t%0.2f')]);
    disp(['Chondritic input H2O is extra     ' num2str(100*dWtH2O,'\t%0.1f') ' % of ' Planet.name '''s mass'])
    disp('')
    if Planet.FeCore
        disp(['wt Pct of the core:    ' num2str(100*Wtcore,'\t%0.2f')]);
        disp(['Chondritic input core matter is extra ' num2str(100*dWtcore,'\t%0.1f') ' % of ' Planet.name '''s mass'])
        disp('')
        disp(['wt Pct S as FeS in core:' num2str(100*XS,'\t%0.2f')]);
        disp(['wt Pct core S in Europa:' num2str(100*WtS,'\t%0.2f')]);
        disp(['Chondritic input S is extra      ' num2str(100*dWtS,'\t%0.1f') ' % of ' Planet.name '''s mass'])
        disp('')
    end
    end
    if Planet.FeCore
    disp(['R Fe (km):             ' num2str(R_Fe_mean_m*1e-3,'\t%0.0f')])
    disp(['span R Fe (km):        ' num2str(R_Fe_range_m*1e-3,'\t%0.0f')])
    disp(['R sil (km):            ' num2str(R_sil_mean_m*1e-3,'\t%0.0f')])
    disp(['span R sil (km):       ' num2str(R_sil_range_m*1e-3,'\t%0.0f')])
    end
end


%% Iterative Calculation for Adding Details to the Seismic Profile
% consider, if appropriate, a convective ice shell of the same thickness
% this introduces an error in the gravity profile and thus the moment of inertia, since the overlying mass will be
% less
%for iT = 1:nTbs
    % Ice I from ConvectionDeschampsSotin above
    
    %was moved up earlier in code
%     if ~isfield(Params,'NO_ICEI_CONVECTION') || Params.NO_ICEI_CONVECTION == false
%         if CONVECTION_FLAG_I(iT) || (isfield(Params,'FORCE_ICEI_CONVECTION') && Params.FORCE_ICEI_CONVECTION == true)
%             %conductive upper layer
%             nConvectIce=n_iceI-nIceIIILithosphere-1;
%             inds_eTBL = find(z_m(iT,1:nConvectIce)<=eTBL_m(iT));
%             Pterm = P_MPa(iT,inds_eTBL)./P_MPa(iT,inds_eTBL(end));
%             T_K(iT,inds_eTBL) = (Tc(iT).^(Pterm)).*(Planet.Tsurf_K.^(1-Pterm));
%             %convective region
%             for iconv = inds_eTBL(end)+1:nConvectIce
%     %             rho = 1000./getVspChoukroun2010(P_MPa(iT,iconv),T_K(iT,iconv-1),2);
%                 rho = getRhoIce(P_MPa(iT,iconv),T_K(iT,iconv-1),1);
%                 try
%                     [Cp,alpha_K] = getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),1);
%                 catch
%                     warning('Seafreeze couldn''t get Cp, alpha for ice. Is it installed?');
%                 end
%     %             aK = 1.56e-4;
%                 T_K(iT,iconv) = T_K(iT,iconv-1)+alpha_K*T_K(iT,iconv)/Cp/rho*deltaP*1e6;
%             end
%            % conductive lower layer
%            inds_deltaTBL = find(z_m(iT,1:nConvectIce)>=z_m(iT,nConvectIce)-deltaTBL_m(iT));
%            T_K(iT,inds_deltaTBL) = (Planet.Tb_K(iT).^(P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT))).*(T_K(iT,inds_deltaTBL(1)-1).^(1-P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT)));
%            rho_kgm3(iT,inds_deltaTBL) = getRhoIce(P_MPa(iT,inds_deltaTBL),T_K(iT,inds_deltaTBL),1);
%
%            if find(phase(iT,:)>1)
%     %         indVI = find(phase(iT,:)==6);
%     %         Ttop = T_K(iT,indsVI(iT));
%     %         Tbottom = zVI_m
%     %         [Q_Wm2(iT),deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),CONVECTION_FLAG_I(iT)]=...
%     %         ConvectionHPIceKalousova2018(Ttop,,PbI_MPa(iT)/2,zb_outerIce_m(iT),g_ms2(iT,1),2);
%            end
%
%         %else % We take care of this case well above in order to be able
%               % to exit sooner if cfg.SKIP_PROFILES=1
%         end
%     else
%         Q_Wm2(iT) = Planet.Qmantle_Wm2(1);
%     end
    
%end
%% Calculate Ice velocities
if Params.CALC_NEW_SOUNDSPEEDS
    disp('Computing sound speeds')
%     velsIce = iceVelsGagnon1990(P_MPa,T_K);
    [vfluid_kms,Ksfluid_GPa] = deal(zeros(nTbs,nsteps));
    for iT = 1:nTbs
        tic
        out=Helgerud_sI(P_MPa(iT,:)', T_K(iT,:)');
        velsIce.Vclathl_kms(iT,:) = 1e-3*out.Vp;
        velsIce.Vclatht_kms(iT,:) = 1e-3*out.Vs;
        velsIce.Ksclath_GPa(iT,:) = 1e-3*out.K;
        velsIce.Gsclath_GPa(iT,:) = 1e-3*out.shear;
        
        out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'Ih');
        velsIce.VIl_kms(iT,:) = 1e-3*out.Vp;
        velsIce.VIt_kms(iT,:) = 1e-3*out.Vs;
        velsIce.KsI_GPa(iT,:) = 1e-3*out.Ks;
        velsIce.GsI_GPa(iT,:) = 1e-3*out.shear;

        out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'III');
        velsIce.VIIIl_kms(iT,:) = 1e-3*out.Vp;
        velsIce.VIIIt_kms(iT,:) = 1e-3*out.Vs;
        velsIce.KsIII_GPa(iT,:) = 1e-3*out.Ks;
        velsIce.GsIII_GPa(iT,:) = 1e-3*out.shear;

        out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'V');
        velsIce.VVl_kms(iT,:) = 1e-3*out.Vp;
        velsIce.VVt_kms(iT,:) = 1e-3*out.Vs;
        velsIce.KsV_GPa(iT,:) = 1e-3*out.Ks;
        velsIce.GsV_GPa(iT,:) = 1e-3*out.shear;
        
        out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'VI');
        velsIce.VVIl_kms(iT,:) = 1e-3*out.Vp;
        velsIce.VVIt_kms(iT,:) = 1e-3*out.Vs;
        velsIce.KsVI_GPa(iT,:) = 1e-3*out.Ks;
        velsIce.GsVI_GPa(iT,:) = 1e-3*out.shear;
        
        
        if POROUS_ICE
            % only affects top ice/clath layers
            % correction for porosity
            por_in(iT).p=P_MPa(iT,1:(n_clath(iT)+n_iceI(iT)))*1e-3;
            por_in(iT).t = T_K(iT,1:(n_clath(iT)+n_iceI(iT)));
            por_in(iT).den = rho_kgm3(iT,1:(n_clath(iT)+n_iceI(iT)));
            try
                por_in(iT).phi_surface=Planet.phi_surface;
            end
            
            por_in(iT).vp = velsIce.Vclathl_kms(iT,1:(n_clath(iT)+n_iceI(iT)));
            por_in(iT).vs = velsIce.Vclatht_kms(iT,1:(n_clath(iT)+n_iceI(iT)));
            por_vel_out(iT) = porosity_correction_ice(por_in(iT),por(iT,:));
            velsIce.Vclathl_kms(iT,1:(n_clath(iT)+n_iceI(iT))) = por_vel_out(iT).vp;
            velsIce.Vclatht_kms(iT,1:(n_clath(iT)+n_iceI(iT))) = por_vel_out(iT).vs;
            
            
            por_in(iT).vp = velsIce.VIl_kms(iT,1:(n_clath(iT)+n_iceI(iT)));
            por_in(iT).vs = velsIce.VIt_kms(iT,1:(n_clath(iT)+n_iceI(iT)));
            por_vel_out(iT) = porosity_correction_ice(por_in(iT),por(iT,:));
            velsIce.VIl_kms(iT,1:(n_clath(iT)+n_iceI(iT))) = por_vel_out(iT).vp;
            velsIce.VIt_kms(iT,1:(n_clath(iT)+n_iceI(iT))) = por_vel_out(iT).vs;
        end
    
        ir = find(phase(iT,:)==0);
        vfluid_kms(iT,ir) = fluidSoundSpeeds(P_MPa(iT,ir),T_K(iT,ir),wo,Planet.Ocean.comp);
        ireal = find(~isnan(vfluid_kms(iT,:)));
        if length(ireal)<length(ir)
            disp('warning: extrapolating fluid sound speeds; this is seawater above 120 MPa, right?')
            vfluid_kms(iT,ir) = interp1(P_MPa(iT,ireal),vfluid_kms(iT,ireal),P_MPa(iT,ir),'linear','extrap');
        end
%         for ir = find(phase(iT,:)==0)
%             vfluid_kms(iT,ir) = fluidSoundSpeeds(P_MPa(iT,ir),T_K(iT,ir),wo,Planet.Ocean.comp);
%         end
        toc
        Ksfluid_GPa(iT,ir)=1e-3*rho_kgm3(iT,ir).*vfluid_kms(iT,ir).^2;
        vfluid_kms(iT,vfluid_kms(iT,:)==0)=NaN;
        Ksfluid_GPa(iT,vfluid_kms(iT,:)==0)=NaN;
        velsIce.VIl_kms(iT,phase(iT,:)~=1) =NaN;
        velsIce.VIt_kms(iT,phase(iT,:)~=1) =NaN;
        velsIce.VIIl_kms(iT,phase(iT,:)~=2) =NaN;
        velsIce.VIIt_kms(iT,phase(iT,:)~=2) =NaN;
        velsIce.VIIIl_kms(iT,phase(iT,:)~=3) =NaN;
        velsIce.VIIIt_kms(iT,phase(iT,:)~=3) =NaN;
        velsIce.VVl_kms(iT,phase(iT,:)~=5) =NaN;
        velsIce.VVt_kms(iT,phase(iT,:)~=5) =NaN;
        velsIce.VVIl_kms(iT,phase(iT,:)~=6) =NaN;
        velsIce.VVIt_kms(iT,phase(iT,:)~=6) =NaN;
        velsIce.Vclathl_kms(iT,phase(iT,:)~=30) =NaN;
        velsIce.Vclatht_kms(iT,phase(iT,:)~=30) =NaN;
    end
    save(fullfile([datpath savefile '_Vels' '_pp' num2str(iT)]),'Ksfluid_GPa','velsIce','vfluid_kms');
else
    try	
        load(fullfile([datpath savefile '_Vels' '_pp' num2str(iT)]),'Ksfluid_GPa','velsIce','vfluid_kms');	
    catch	
        error(['ERROR: A soundspeeds file was not found for ' char(Planet.salt) '. Re-run with CALC_NEW_SOUND=1 to generate the needed file.'])	
    end
end

if isfield(Seismic,'mantleEOS')
    mantle = loadMantleEOS(Seismic.mantleEOS);
end
% mprops = load(Seismic.mantleEOS);
if isfield(Seismic,'mantleEOS_dry')
    mantleDry = loadMantleEOS(Seismic.mantleEOS_dry);
end
    

for iT = 1:nTbs
    mean_rho_ocean(iT) = mean(rho_kgm3(iT,n_iceI(iT)+1:indSil(iT)));
end
    disp(['rho ocean (kg/m3) ' num2str(mean_rho_ocean,'\t%0.0f')])


if isfield(Seismic,'SMOOTH_ROCK') && Params.SMOOTH_VROCK
    ncr = Seismic.SMOOTH_VROCK;
    mantle.vp = smooth2a(mantle.vp,ncr,ncr);
    mantle.vs = smooth2a(mantle.vs,ncr,ncr);
end

%% Core
% Iron properties
    %gamma (fcc) Fe in Table 3 of C2006. derivatives are  rounded off values for alpha (bcc) Fe
    % gamma Fe should be the correct phase between 1394 C (iron allotropes
    % wiki; 1667 K) and 912 C (1185 K).  Minimum Tcore in C2006 is
    % 1200K.
%     rhoo_iron_kgm3 = 8e3;
if Planet.FeCore
    rhoo_iron_kgm3 = rho_Fe_kgm3; %determined above
else
    rhoo_iron_kgm3 = nan;
end
    if Planet.xFeS<0.05 %Weight fraction of sulfur in the core
        alpha_iron_K = 5e-5;
    %     Bulk modulus and derivatives?
        Kso_iron_GPa = 156;
        dKsdP_iron = 5;
        dKsdT_PaK = -0.040;
    %     shear modulus and derivatives
        Go_iron_GPa = 76.5;
        dGdP_iron = 2;
        dGdT_iron_PaK = -0.023;
    else
        rhoo_iron_kgm3 = 5150-50*(100*Planet.xFeS-10); %Only really valid for Planet.xFeS<0.2 (Cammarano et al., 2006, Table 3)
        alpha_iron_K = 9.2e-5;
    %     Bulk modulus and derivatives?
        Kso_iron_GPa = 53.2-2*(100*Planet.xFeS-10);
        dKsdP_iron = 4.66;
        dKsdT_PaK = 0;
    %     shear modulus and derivatives
        Go_iron_GPa = 0;
        dGdP_iron = 0;
        dGdT_iron_PaK = 0;
    end

[Pcore_MPa,Tcore_K,rhocore_kgm3,M_above_core,VP_core_kms,VS_core_kms,g_ms2_core,Ks_iron_GPa,G_iron_GPa] = deal(zeros(nTbs,Params.nsteps_core));
clear g_ms2_core M_above_core
    
for iT = nTbs:-1:1  % Do this loop in decreasing order to avoid preallocating por_in, por_out structs
    [Pmantle_MPa,rho_mantle_kgm3,M_above_mantle,VP_mantle_kms,VS_mantle_kms,g_ms2_sil] = deal(zeros(1,nsteps_mantle(iT)));
    Tmantle_K = T_K(iT,indSil(iT));
    r_mantle_m = linspace(R_sil_mean_m(iT),R_Fe_mean_m(iT),nsteps_mantle(iT));
    
    g_ms2_sil = g_ms2(iT,C2mean(iT))*ones(1,nsteps_mantle(iT));
    MantleHeat = Planet.Qmantle_Wm2(iT)*4*pi*Planet.R_m^2+Planet.QHmantle;
    rhom_rough = 3000;
    alpha_rough = 0.2e-4;
    Cp_rough = 2e6;
    Kappa_rough = 1e-6;
    nu_mantle = 1e21; % mantle viscosity in Pa S, a common number for Earth?s mantle
    DeltaT = 800; % temperature differential in K
    
    % cold case
    if Planet.FeCore
        Tmantle_K = conductiveMantleTemperature(r_mantle_m,R_Fe_mean_m(iT),R_sil_mean_m(iT),Planet.kr_mantle,Planet.rho_sil_withcore_kgm3,Tmantle_K,MantleHeat,0);
        MantleMass = 4/3*pi*(R_sil_mean_m(iT)^3-R_Fe_mean_m(iT)^3);
        Dm = R_sil_mean_m(iT)-R_Fe_mean_m(iT);
    else
        Tmantle_K = conductiveMantleTemperature(r_mantle_m,0,R_sil_mean_m(iT),Planet.kr_mantle,rho_sil_kgm3(iT,C2mean(iT)),Tmantle_K,MantleHeat,0);
        MantleMass = 4/3*pi*(R_sil_mean_m(iT)^3);
        Dm = R_sil_mean_m(iT);
    end    
    
    Ra_mantle = g_ms2_sil(1)*rhom_rough*alpha_rough*MantleHeat/MantleMass*Dm^5/nu_mantle/Kappa_rough^2/Cp_rough; % bunge 1997
    % parmentier and sotin 1994
    Ra_mantle = alpha_rough*rhom_rough*g_ms2_sil(1)*DeltaT*Dm^3/Kappa_rough/nu_mantle;% as per Hussmann and Spohn 2004

    % COMMON FUNCTIONALITY HERE IS TO PROPAGATE PRESSURE, GRAVITY AND
    % DENSITY DOWNWARD FOR A CONDUCTIVELY COOLING LAYER WITH NO SOLID STATE
    % CONVECTION
%     g_ms2_sil(1) = g_ms2(iT,C2mean(iT));
    M_above_mantle(1) = M_above_kg(iT,C2mean(iT));
    Pmantle_MPa(1) = P_MPa(iT,C2mean(iT));
    rhofn_kgm3 = mantle.rho_fn;
    VPfn_kms = mantle.vp_fn;
    VSfn_kms = mantle.vs_fn;
    rho_mantle_kgm3(1) = rhofn_kgm3(Pmantle_MPa(1),Tmantle_K(1));
    VP_mantle_kms(1) = VPfn_kms(Pmantle_MPa(1),Tmantle_K(1));
    VS_mantle_kms(1) = VSfn_kms(Pmantle_MPa(1),Tmantle_K(1));

    for ij = 2:nsteps_mantle(iT)
%         if r_mantle_m(ij-1)>0.3*Planet.R_m
%             g_ms2_sil(ij)=(Gg*(Planet.M_kg-M_above_mantle(ij-1))/r_mantle_m(ij-1)^2);
%         else
%             g_ms2_sil(ij) = g_ms2_sil(ij-1);
%         end
        Pmantle_MPa(ij) = Pmantle_MPa(ij-1)+1e-6*(rho_mantle_kgm3(ij-1))*g_ms2_sil(ij)*(r_mantle_m(ij-1)-r_mantle_m(ij));
        M_above_mantle(ij) = M_above_mantle(ij-1)+4/3*pi*(r_mantle_m(ij-1)^3-r_mantle_m(ij)^3)*rho_mantle_kgm3(ij-1);
        
        if isfield(Seismic,'mantleEOS_dry') && rho_mantle_kgm3(ij-1)>rhofn_kgm3(Pmantle_MPa(ij),Tmantle_K(ij))
            rhofn_kgm3 = mantleDry.rho_fn;
            VPfn_kms = mantleDry.vp_fn;
            VSfn_kms = mantleDry.vs_fn;
        end
         rho_mantle_kgm3(ij) = rhofn_kgm3(Pmantle_MPa(ij),Tmantle_K(ij));
         VP_mantle_kms(ij) = VPfn_kms(Pmantle_MPa(ij),Tmantle_K(ij));
         VS_mantle_kms(ij) = VSfn_kms(Pmantle_MPa(ij),Tmantle_K(ij));
         if VP_mantle_kms(ij)<0.1 || isnan(VP_mantle_kms(ij))
             VP_mantle_kms(ij) = NaN;
         end
         if VS_mantle_kms(ij)<0.1 || isnan(VS_mantle_kms(ij))
             VS_mantle_kms(ij) = NaN;
         end
         if Pmantle_MPa(ij)<0.1 || isnan(Pmantle_MPa(ij))
             Pmantle_MPa(ij) = NaN;
         end

    end
    
    if POROUS
    % correction for porosity
        por_in(iT).p=Pmantle_MPa*1e-3;
        por_in(iT).t = Tmantle_K;
        por_in(iT).den = rho_mantle_kgm3;
        por_in(iT).vp = VP_mantle_kms;
        por_in(iT).vs = VS_mantle_kms;
        if isfield(Planet,'phi_surface')
            por_out(iT) = porosity_correction(por_in(iT),Planet.phi_surface);
        else
            por_out(iT) = porosity_correction(por_in(iT));
        end
        permeability = por_out(iT).per;
        rho_mantle_kgm3 = por_out(iT).den;
        VP_mantle_kms = por_out(iT).vp;
        VS_mantle_kms = por_out(iT).vs;
        disp(['Average Porosity: ' num2str(mean(por_out(iT).por))])
        disp(['Porosity: ' num2str(por_out(iT).por)])
        %recalculate mass above in light of reduced density. This should
        %really be done in a recursive way that acknowledges the reduced
        %overburden pressure.
        for ij = 2:nsteps_mantle(iT)
%         if r_mantle_m(ij-1)>0.3*Planet.R_m
%             g_ms2_sil(ij)=(Gg*(Planet.M_kg-M_above_mantle(ij-1))/r_mantle_m(ij-1)^2);
%         else
%             g_ms2_sil(ij) = g_ms2_sil(ij-1);
%         end
%         Pmantle_MPa(ij) = Pmantle_MPa(ij-1)+1e-6*(rho_mantle_kgm3(ij-1))*g_ms2_sil(ij)*(r_mantle_m(ij-1)-r_mantle_m(ij));
          M_above_mantle(ij) = M_above_mantle(ij-1)+4/3*pi*(r_mantle_m(ij-1)^3-r_mantle_m(ij)^3)*rho_mantle_kgm3(ij-1);
        end
    end
    
    mtest = find(M_above_mantle>Planet.M_kg);
    if mtest
        disp(['Exceeded planet mass at mantle radius of ' num2str(1e-3*r_mantle_m(mtest(1)), '%0.0f'), ' km']);
        
    end
    
    if find(g_ms2_sil<0)
        disp('too-large input density (probably in the mantle) resulted in negative gravity; computing pressure and gravity from density instead')
    end
    
%%  insert a low velocity layer in upper part of mantle; this was useful to address some quirks of AxiSEM
%     indsLow = find(r_mantle_m(1)-r_mantle_m<=5000);strLow = 'LowVUpperMantle5km';
%     indsLow = find(r_mantle_m(1)-r_mantle_m<=20000); strLow = 'LowVUpperMantle20km';
% if Seismic.LOWV_MANTLE
%     indsLow = find(r_mantle_m(1)-r_mantle_m<=30000); strLow = 'LowVUpper30kmMantle';
%     VP_mantle_kms(indsLow) = 4/7* VP_mantle_kms(indsLow);
%     VS_mantle_kms(indsLow) = 4/7* VS_mantle_kms(indsLow);
% end
strLow = '';
%%
    % AS ABOVE, NOW FOR THE METALLIC CORE; THESE VALUES AREN'T USED IF PLANET.FECORE IS FALSE
%     g_ms2_core(iT,1) = g_ms2_sil(nsteps_mantle(iT)); % depth dependent
%     gravity was used previously, but this can result in negative gravity
%     approaching the center of the core due to imperfect matching of the
%     moment of inertia.
%     M_above_core(iT,1) = M_above_mantle(nsteps_mantle(iT));
    if Planet.FeCore
        M_above_core(iT,1) = M_above_mantle(nsteps_mantle(iT));
        g_ms2_core(iT,:)=(Gg*(Planet.M_kg-M_above_core(iT))/r_core_m(iT)^2)*ones(1,Params.nsteps_core);

        Tcore_K(iT,:) = linspace(Tmantle_K(end),1.01*Tmantle_K(end),Params.nsteps_core);
        Pcore_MPa(iT,1) = Pmantle_MPa(nsteps_mantle(iT));
        if isfield(Seismic,'coreEOS')
            core = loadMantleEOS(Seismic.coreEOS);
            rhocore_kgm3(iT,1) = core.rho_fn(Pcore_MPa(iT,1),Tcore_K(iT,1));
            
            for ij = 2:Params.nsteps_core
        %          g_ms2_core(iT,ij)=(Gg*(Planet.M_kg-M_above_core(iT,ij-1))/r_core_m(iT,ij-1)^2);
        %          Pcore_MPa(iT,ij) = Pcore_MPa(iT,ij-1)+1e-6*rhocore_kgm3(iT,ij-1)*g_ms2_core(iT,ij)*(r_core_m(iT,ij-1)-r_core_m(iT,ij));
                 Pcore_MPa(iT,ij) = Pcore_MPa(iT,ij-1)+1e-6*rhocore_kgm3(iT,ij-1)*g_ms2_core(iT)*(r_core_m(iT,ij-1)-r_core_m(iT,ij));
                 rhocore_kgm3(iT,ij) = core.rho_fn(Pcore_MPa(iT,ij),Tcore_K(iT,ij));
    %              rhocore_kgm3(iT,ij) = rhocore_kgm3(iT,ij-1)*(1+1./(1e-3*Pcore_MPa(iT,ij)*Ks_iron_GPa(iT,ij)));
        %          M_above_core(iT,ij) = M_above_core(iT,ij-1)+4/3*pi*(r_core_m(iT,ij-1)^3-r_core_m(iT,ij)^3)*rhocore_kgm3(iT,ij-1);
            end
            Ks_iron_GPa(iT,:) = bar2GPa*core.Ks_fn(Pcore_MPa(iT,:),Tcore_K(iT,:));
            G_iron_GPa(iT,:) = bar2GPa*core.Gs_fn(Pcore_MPa(iT,:),Tcore_K(iT,:));
            VS_core_kms(iT,:) = core.vs_fn(Pcore_MPa(iT,:),Tcore_K(iT,:));
            VP_core_kms(iT,:) = core.vp_fn(Pcore_MPa(iT,:),Tcore_K(iT,:));
        else
            Ks_iron_GPa(iT,1) = Kso_iron_GPa+1e-3*Pcore_MPa(iT,1)*dKsdP_iron+1e-9*Tcore_K(iT,1)*dKsdT_PaK;
            G_iron_GPa(iT,1) = Go_iron_GPa+1e-3*Pcore_MPa(iT,1)*dGdP_iron+1e-9*Tcore_K(iT,1)*dGdT_iron_PaK;
            rhocore_kgm3(iT,:) = rhoo_iron_kgm3;

            for ij = 2:Params.nsteps_core
        %          g_ms2_core(iT,ij)=(Gg*(Planet.M_kg-M_above_core(iT,ij-1))/r_core_m(iT,ij-1)^2);
        %          Pcore_MPa(iT,ij) = Pcore_MPa(iT,ij-1)+1e-6*rhocore_kgm3(iT,ij-1)*g_ms2_core(iT,ij)*(r_core_m(iT,ij-1)-r_core_m(iT,ij));
                 Pcore_MPa(iT,ij) = Pcore_MPa(iT,ij-1)+1e-6*rhocore_kgm3(iT,ij-1)*g_ms2_core(iT)*(r_core_m(iT,ij-1)-r_core_m(iT,ij));
                 Ks_iron_GPa(iT,ij) = Kso_iron_GPa+1e-3*Pcore_MPa(iT,ij)*dKsdP_iron+1e-9*Tcore_K(iT,ij)*dKsdT_PaK;
                 G_iron_GPa(iT,ij) = Go_iron_GPa+1e-3*Pcore_MPa(iT,ij)*dGdP_iron+1e-9*Tcore_K(iT,ij)*dGdT_iron_PaK;
    %              rhocore_kgm3(iT,ij) = rhocore_kgm3(iT,ij-1)*(1+1./(1e-3*Pcore_MPa(iT,ij)*Ks_iron_GPa(iT,ij)));
        %          M_above_core(iT,ij) = M_above_core(iT,ij-1)+4/3*pi*(r_core_m(iT,ij-1)^3-r_core_m(iT,ij)^3)*rhocore_kgm3(iT,ij-1);
            end
            VS_core_kms(iT,:) = 1e-3*sqrt(G_iron_GPa(iT,:)*1e9./rhocore_kgm3(iT,:));
            VP_core_kms(iT,:) = 1e-3*sqrt(Ks_iron_GPa(iT,:)*1e9./rhocore_kgm3(iT,:)+4/3*(VS_core_kms(iT,:)*1e3).^2);
        end
    end

%     kluge to fix NaN's in the mantle densities and sound speeds
    rho_mantle_kgm3(isnan(rho_mantle_kgm3))=0;
    VS_mantle_kms(isnan(VS_mantle_kms))=0;
    VP_mantle_kms(isnan(VP_mantle_kms))=0;
    
    
    interior(iT).g_ms2 = g_ms2_sil;
    interior(iT).r_mantle_m = r_mantle_m;
    interior(iT).rho_mantle_kgm3 = rho_mantle_kgm3;
    interior(iT).Tmantle_K = Tmantle_K;
    interior(iT).Pmantle_MPa = Pmantle_MPa;
    interior(iT).VS_mantle_kms = interp1nan(Pmantle_MPa,VS_mantle_kms);
    interior(iT).VP_mantle_kms = interp1nan(Pmantle_MPa,VP_mantle_kms);
    if POROUS
        interior(iT).permeability = permeability;
    end
    if isfield(Seismic,'mfluids') % still under development. A goal is to more realistically track the loss of fluids along the geotherm.
        mfluids = loadMantleFluids(Seismic.mfluids);
        mphases = loadMantlePhases(Seismic.mphases);
        mvolumes = loadMantlePhases(Seismic.mvolumes);

        mfnames = fieldnames(mfluids);
        mpnames = fieldnames(mphases);
        mvnames = fieldnames(mvolumes);
        for im = 1:length(mfnames) % count only the names, not the interpolating functions or p or t
            if ~endsWith(mfnames{im},'_fn') && ~(strcmp(mfnames{im},'p') || strcmp(mfnames{im},'t'))
                interior(iT).(['fluid_' mfnames{im}]) = mfluids.(mfnames{im+1})(Pmantle_MPa,Tmantle_K);
            end
        end
        for im = 1:length(mpnames) % count only the names, not the interpolating functions or p or t
            if ~endsWith(mpnames{im},'_fn') && ~(strcmp(mpnames{im},'p') || strcmp(mpnames{im},'t'))
                interior(iT).(['wt_' mpnames{im}]) = mphases.(mpnames{im+1})(Pmantle_MPa,Tmantle_K); %wt pct of each constituent
            end
        end
        masscheck = 0;
        for im = 1:length(mvnames) % count only the names, not the interpolating functions or p or t
            if ~endsWith(mvnames{im},'_fn') && ~(strcmp(mvnames{im},'p') || strcmp(mvnames{im},'t'))
                disp(mvnames{im})
                interior(iT).(['vol_' mvnames{im}]) = mvolumes.(mvnames{im+1})(Pmantle_MPa,Tmantle_K); %vol pct of each constituent
                interior(iT).(['rho_' mvnames{im}]) = interior(iT).(['wt_' mpnames{im}])./interior(iT).(['vol_' mvnames{im}]).*rho_mantle_kgm3; %dens fraction of each constituent
                interior(iT).(['mass_' mvnames{im}]) = 4/3*pi*interior(iT).(['rho_' mvnames{im}]).*-gradient(r_mantle_m.^3); %mass in each spherical layer of the mantle
            end
        end
    end
    if isfield(mantle,'Ks_fn')
    interior(iT).Ks_GPa = bar2GPa*mantle.Ks_fn(Pmantle_MPa,Tmantle_K);
    interior(iT).Ks_GPa = interp1nan(Pmantle_MPa,interior(iT).Ks_GPa);
    interior(iT).Gs_GPa = bar2GPa*mantle.Gs_fn(Pmantle_MPa,Tmantle_K);
    interior(iT).Gs_GPa = interp1nan(Pmantle_MPa,interior(iT).Gs_GPa);
    end
        % compute seismic attenuation, as per C2006,
    % Rb = 8.314462175; % J/K/mol
    Tm_K = 273.15+Tm_p_Hirschmann2000(1e-3*Pmantle_MPa); %input in GPa
    Hp_mantle = Seismic.g_aniso_mantle*Tm_K;
    interior(iT).QS_overfgamma = Seismic.B_aniso_mantle*exp(Seismic.gamma_aniso_mantle*Hp_mantle./Tmantle_K);
%%
    if POROUS && ~cfg.SKIP_PROFILES
        if Params.HOLD
            set(0, 'CurrentFigure', figs.porP);
            set(gcf, 'Name', lbl.porvP)
        else
            set(0, 'CurrentFigure', figs.porP(iT));
            set(gcf, 'Name', [lbl.porvP ' Tb = ' num2str(Planet.Tb_K(iT))])
            clf;
        end
        hold on
%         plot(Pmantle_MPa,[por_in(iT).vp]-[por_out(iT).vp])
%         plot(Pmantle_MPa,por_in(iT).vs-por_out(iT).vs,'--')
        plot(Pmantle_MPa,por_out(iT).por*100,'LineWidth',cfg.LW_std)
%         plot(Pmantle_MPa,(por_in(iT).den-por_out(iT).den)./por_in(iT).den*100)
%         plot(Pmantle_MPa,(por_in(iT).vp-por_out(iT).vp)./por_in(iT).vp*100)
%         plot(Pmantle_MPa,(por_in(iT).vs-por_out(iT).vs)./por_in(iT).vs*100,'--')
        xlabel([math 'P_{' nm 'rock}' nm ' (MPa)'],'Fontsize',lbl.mltext)
%         ylabel('v-v_{porous} (m s^{-1})');
%         ylabel('$\frac{X-X_{porous}}{X}$ (\%)');
        ylabel(['Porosity ' math '\phi' nm ' (%)'],'Fontsize',lbl.mltext);
%         legend('\phi','\rho','V_P','V_S')
        box on
        axis tight
        
        if Params.HOLD
            set(0, 'CurrentFigure', figs.porR);
            set(gcf, 'Name', lbl.porvR)
        else
            set(0, 'CurrentFigure', figs.porR(iT));
            set(gcf, 'Name', [lbl.porvR ' Tb = ' num2str(Planet.Tb_K(iT))])
            clf;
        end
        hold on
%         plot(Pmantle_MPa,[por_in(iT).vp]-[por_out(iT).vp])
%         plot(Pmantle_MPa,por_in(iT).vs-por_out(iT).vs,'--')
        hl = plot(por_out(iT).por*100,r_mantle_m*1e-3,'LineWidth',cfg.LW_std);

%         plot((por_in(iT).den-por_out(iT).den)./por_in(iT).den*100,r_mantle_m*1e-3)
%         plot((por_in(iT).vp-por_out(iT).vp)./por_in(iT).vp*100,r_mantle_m*1e-3)
%         plot((por_in(iT).vs-por_out(iT).vs)./por_in(iT).vs*100,r_mantle_m*1e-3,'--')
%         set(gca,'ydir','reverse');
        ylabel([math 'r_{' nm 'rock}' nm ' (km)'],'Fontsize',lbl.mltext)
%         ylabel('v-v_{porous} (m s^{-1})');
%         xlabel('$\frac{X-X_{porous}}{X}$ (\%)');
        xlabel(['Porosity ' math '\phi' nm ' (%)'],'Fontsize',lbl.mltext);

%         legend('\phi','\rho','V_P','V_S')
        box on
        axis tight
        
        if Params.HOLD
            set(0, 'CurrentFigure', figs.perm);
            set(gcf, 'Name', lbl.perme)
        else
            set(0, 'CurrentFigure', figs.perm(iT));
            set(gcf, 'Name', [lbl.perme ' Tb = ' num2str(Planet.Tb_K(iT))])
            clf;
        end
        hold on
%         plot(Pmantle_MPa,[por_in(iT).vp]-[por_out(iT).vp])
%         plot(Pmantle_MPa,por_in(iT).vs-por_out(iT).vs,'--')
%         plot(squeeze(log10(permeability)),r_mantle_m'*ones(1,length(permeability(1,1,:)))*1e-3)
        hl1 = plot(squeeze(log10(permeability(1,:,1))),r_mantle_m*1e-3,'LineWidth',cfg.LW_std);
        hl2 = plot(squeeze(log10(permeability(1,:,2))),r_mantle_m*1e-3,'--','LineWidth',cfg.LW_std);
        hl3 = plot(squeeze(log10(permeability(1,:,3))),r_mantle_m*1e-3,'LineWidth',cfg.LW_std);
        hl4 = plot(squeeze(log10(permeability(1,:,4))),r_mantle_m*1e-3,'LineWidth',cfg.LW_std);
        hl5 = plot(squeeze(log10(permeability(1,:,5))),r_mantle_m*1e-3,'LineWidth',cfg.LW_std);

        legend({'Crust in general','Upper crust in general','Low permeability upper crust','Disturbed crust','Oceanic crust'},'Fontsize',lbl.smtext);

%         set(gca,'ydir','reverse');
        ylabel([math 'r_{' nm 'mantle}' nm ' (km)'],'Fontsize',lbl.mltext)
%         ylabel('v-v_{porous} (m s^{-1})');
        xlabel([math '\rm log_{10}' nm ' permeability'],'Fontsize',lbl.mltext);
        box on
        
        % Save porosity figures later because we need more info to create 'thissavestr'
    end % POROUS && ~cfg.SKIP_PROFILES
end

if cfg.DISP_LAYERS
    %% create the table of depths and heat flux values
    % multirow{4}{*}{5 wt\%} &Ice Ih &-&141& 127  & 96 & 63  & 24 \\
    % & Liquid &-&63&92&128&244&539\\
    % & Ice III &-&29&-&-&-&-\\
    % & Ice V &-&158&134&49&-&-\\
    % & Ice VI &-&409&411&411&355&237\\

    Tb_str = getTableStr(Planet.Tb_K,Planet.Tb_K);
    Qb_str = getTableStr(Planet.Tb_K,round(Qb*1e3));
    Qc_str = getTableStr(Planet.Tb_K,round(Q_Wm2*1e3));
    nu_str = getTableStr(Planet.Tb_K,log10(nu));
    dzI_str = getTableStr(Planet.Tb_K,round(zI_m*1e-3));
    dz_ocean_m_str = getTableStr(Planet.Tb_K,round(dz_ocean_m_m*1e-3));
    dzIII_str = getTableStr(Planet.Tb_K,round(dzIII_m*1e-3));
    dzV_str = getTableStr(Planet.Tb_K,round(dzV_m*1e-3));
    dzVI_str = getTableStr(Planet.Tb_K,round(dzVI_m*1e-3));
    R_Fe_str = getTableStr(Planet.Tb_K,round(R_Fe_mean_m*1e-3));
    R_sil_str = getTableStr(Planet.Tb_K,round((R_sil_mean_m)*1e-3));

    disp('\hline')
    % for latex environments supporting multirow (not AGU)
    % if wo==0
    %     disp(['\multirow{10}{*}{\begin{tabular}{c}Water\end{tabular}}'])
    % else
    %     switch Planet.Ocean.comp
    %         case 'MgSO4'
    %             disp('\multirow{10}{*}{\begin{tabular}{c} MgSO$_{4}$ \\' );
    %             disp([num2str(wo) 'wt\% \end{tabular}}']);
    %         case 'NH3'
    %             disp('\multirow{10}{*}{\begin{tabular}{c} NH$_{3}$\\')
    %             disp([num2str(wo) 'wt\% \end{tabular}}']);
    %         case 'Seawater'
    %             disp('\multirow{10}{*}{\begin{tabular}{c} Seawater\\')
    %             disp([num2str(wo) 'g/kg \end{tabular}}']);
    %     end
    % end
    if wo==0
        concstr{1} = 'Water';
        concstr{2} = '';
    else
        switch Planet.Ocean.comp
            case 'MgSO4'
               concstr{1} = 'MgSO$_{4}$';
               concstr{2}= [num2str(wo) 'wt\%'];
            case 'NH3'
               concstr{1}='NH$_{3}$';
               concstr{2}= [num2str(wo) 'wt\%'];
            case 'Seawater'
               concstr{1}= 'Seawater';
               concstr{2}= [num2str(wo) 'g/kg'];
        end
    end
    
    for iT = 1:nTbs
        dV = -gradient([r_m(iT,1:indSil(iT)) interior(iT).r_mantle_m].^3);
        Planet.Mcomputed_kg(iT) = 4/3*pi*sum(dV.*[rho_kgm3(iT,1:indSil(iT)) interior(iT).rho_mantle_kgm3]);
    end
    % Preallocate
    [rhommodel, rhom, meanphivec] = deal(zeros(1,nTbs));

    disp(['Planet Mass (kg): ' num2str(Planet.M_kg)])
    disp(['Computed Mass (kg): ' num2str(Planet.Mcomputed_kg,' %0.4g')]);
    disp(['Input C/MR2: ' num2str(Planet.Cmeasured)])
    if Planet.FeCore
        for iT = 1:nTbs
            disp(['Computed C/MR2 for Tb=' num2str(Tb_K(iT)) 'K: ' num2str(C2(iT,C2mean(iT))/MR2) ' (neighboring values: ' num2str(C2(iT,C2mean(iT)+1)/MR2) '; ' num2str(C2(iT,C2mean(iT)-1)/MR2) ')'])
            rhommodel(iT) = mean(interior(iT).rho_mantle_kgm3);
        end
        rhorockstr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(Planet.rho_sil_withcore_kgm3) '}'];
        rhorockmstr = getTableStr(Planet.Tb_K,round(rhommodel));
        xFeSstr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(100*Planet.xFeS) '}'];
        rhocorestr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(rho_Fe_kgm3,'%0.0f') '}'];
        if exist('concstr','var')
            disp([concstr{1} '&$\rho_{rock}$ (kg m$^{-3}$)&' rhorockstr '\\']);
            disp([concstr{2} '&$\rho_{rock,model}$ (kg m$^{-3}$)' rhorockmstr{:} '\\']);
        else
            disp(['&$\rho_{rock}$ (kg m$^{-3}$)&' rhorockstr '\\']);
            disp(['&$\rho_{rock,model}$ (kg m$^{-3}$)' rhorockmstr{:} '\\']);
        end
        disp(['&$X_{FeS}$ (\%)&' xFeSstr '\\']);
        disp(['&$\rho_{core}$ (kg m$^{-3}$)&' rhocorestr '\\']);
    else
        disp(['Computed C/MR2 for Tb=' num2str(Tb_K(iT)) 'K: ' num2str(C1(iT,C2mean(iT))/MR2) '  (neighboring values: ' num2str(C1(iT,C2mean(iT)+1)/MR2) '; ' num2str(C1(iT,C2mean(iT)-1)/MR2) ')'])
        for iT = 1:nTbs
            rhom(iT) = mean(rho_sil_kgm3(iT,C2mean(:)));
            rhommodel(iT) = mean(interior(iT).rho_mantle_kgm3);
        end
        rhorockstr = getTableStr(Planet.Tb_K,round(rhom));
        rhorockmstr = getTableStr(Planet.Tb_K,round(rhommodel));
        if exist('concstr','var')
            disp([concstr{1} '&$\rho_{rock}$ (kg m$^{-3}$)' rhorockstr{:} '\\']);
            disp([concstr{2} '&$\rho_{rock,model}$ (kg m$^{-3}$)' rhorockmstr{:} '\\']);
        else
            disp(['&$\rho_{rock}$ (kg m$^{-3}$)&' rhorockstr '\\']);
            disp(['&$\rho_{rock,model}$ (kg m$^{-3}$)' rhorockmstr{:} '\\']);
        end
        disp(['\cline{3-' num2str(nTbs+2) '}'])
    end
    Qradrockstr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(Planet.Qmantle_Wm2*4*pi*Planet.R_m^2/1e9) '}'];
    disp(['&$Q_{rock}$ (GW)&' Qradrockstr '\\']);
    if Planet.QHmantle
        Qtiderockstr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(Planet.QHmantle) '}'];
        disp(['&$Q_{rock,tidal}$ (GW)&' Qtiderockstr '\\']);
    end
    disp(['\cline{3-' num2str(nTbs+2) '}'])
    disp(['&T$_{b}$ (K)  ' Tb_str{:} ' \\'])
    if Planet.POROUS_ROCK
        if isfield(Planet,'phi_surface')
            phisurfstr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(100*Planet.phi_surface) '}'];
            disp(['&$\phi_{1}$   ' phisurfstr ' \\'])
        end
        for iT = 1:nTbs
            meanphivec(iT) = 100*mean(por_out(iT).por);
        end
        meanphistr = getTableStr(Planet.Tb_K,round(meanphivec));
        disp(['&$\bar{\phi}$ (\%)  ' meanphistr ' \\'])
    end
    disp(['&$q_\mathrm{b}$ mW m$^{-2}$  ' Qb_str{:} ' \\'])
    disp(['&$q_\mathrm{c}$ mW m$^{-2}$  ' Qc_str{:} ' \\'])
    disp(['&\log_{10}(\nu$_\mathrm{ice})$ mW m$^{-2}$  ' nu_str{:} ' \\'])
    disp(['&$D_\mathrm{Ih}$ (km)' dzI_str{:} ' \\'])
    disp(['&$D_\mathrm{ocean}$ (km)' dz_ocean_m_str{:} ' \\'])
    disp(['&$D_\mathrm{III}$ (km)' dzIII_str{:} ' \\'])
    disp(['&$D_\mathrm{V}$ (km)' dzV_str{:} ' \\'])
    disp(['&$D_\mathrm{VI}$ (km)' dzVI_str{:} ' \\'])
    disp(['&$R_\mathrm{rock}$ (km)' R_sil_str{:} ' \\'])
    if Planet.FeCore
        R_Fe_range_str = getTableStr(Planet.Tb_K,round(R_Fe_range_m*1e-3));
        R_sil_range_str = getTableStr(Planet.Tb_K,round(R_sil_range_m*1e-3));
        disp(['&R$_\mathrm{core}$ (km)' R_Fe_str{:} ' \\'])
        disp(['&$\Delta$R$_\mathrm{core} (km)$' R_Fe_range_str{:} ' \\'])
        disp(['&$\Delta$R$_\mathrm{mantle}$ (km)' R_sil_range_str{:} ' \\'])
    end
    disp('\hline')
end


%% Construct output for seismic modeling:
        % ice                           ocean                                            mantle         core
% [g_Planet_ms2,P_Planet_MPa,T_Planet_K,r_Planet_m,rho_pPlanet_kgm3,VP_Planet_kms,VS_Planet_kms,QS_overfgamma_Planet,k_S_mPlanet]=...
%     deal(zeros(nTbs,indSil(1)-1+nsteps_mantle(1)+Params.nsteps_core));
%% Plot settings
LineStyle=Params.LineStyle;
ymax = 1.05*Planet.R_m*1e-3;

if Planet.FeCore
    [VP_Planet_kms,VS_Planet_kms,Ks_Planet_GPa,Gs_Planet_GPa,QS_overfgamma_Planet,k_S_mPlanet,phasePlanet] = ...
        deal(nan(nTbs,length([g_ms2(1,1:indSil(1)-1) interior(1).g_ms2 g_ms2_core(1,:)])));
else
    [VP_Planet_kms,VS_Planet_kms,Ks_Planet_GPa,Gs_Planet_GPa,QS_overfgamma_Planet,k_S_mPlanet,phasePlanet] = ...
        deal(nan(nTbs,length([g_ms2(1,1:indSil(1)-1) interior(1).g_ms2])));
end

[g_Planet_ms2, P_Planet_MPa, T_Planet_K, r_Planet_m, rho_pPlanet_kgm3] = deal(zeros(nTbs,nsteps_tot));
[kmean, ktop, ocean_thk, ind_Ih, ind_Obot] = deal(nan(1,nTbs));
for iT = 1:nTbs
    % Grab the indices for the ices, but only for the portion of the
    % interior calculation above the silicate interface
    H2Oinds = 1:indSil(iT)-1;
    indsI = find(phase(iT,H2Oinds)==1);
    indsLiquid = find(phase(iT,H2Oinds)==0);
    indsIII = find(phase(iT,H2Oinds)==3);
    indsV = find(phase(iT,H2Oinds)==5);
    indsVI = find(phase(iT,H2Oinds)==6);
    indsClath = find(phase(iT,H2Oinds)==30);
    
    phasePlanet(iT,indsI) = 1;
    phasePlanet(iT,indsLiquid) = 0;
    phasePlanet(iT,indsIII) = 3;
    phasePlanet(iT,indsV) = 5;
    phasePlanet(iT,indsVI) = 6;
    phasePlanet(iT,indsClath) = 30;
    
    phasePlanet(iT,indSil(iT):indSil(iT)+nsteps_mantle(iT)-1)=50; %reserve 50-99 for different mantle types

    if ~Planet.FeCore
        [thisgcore,thisPcore,thisTcore,thisrcore,thisrhocore,thisVPcore,thisVScore,thisQScore,thisKScore,thisGScore,thiskScore]=deal([]); %placeholders to avoid errors
    else
        thisgcore = g_ms2_core(iT,:);
        thisPcore = Pcore_MPa(iT,:);
        thisTcore = Tcore_K(iT,:);
        thisrcore = r_core_m(iT,:);
        thisrhocore = rhocore_kgm3(iT,:);
        thisVPcore = VP_core_kms(iT,:);
        thisVScore = VS_core_kms(iT,:);
        thisQScore = Seismic.QScore*ones(1,Params.nsteps_core);
        thiskScore = 0*Seismic.QScore*ones(1,Params.nsteps_core);
        thisKScore = Ks_iron_GPa(iT,:);
        thisGScore = G_iron_GPa(iT,:);
        phasePlanet(iT,end-Params.nsteps_core+1:end) = 100;
    end
    g_Planet_ms2(iT,:) = [g_ms2(iT,1:indSil(iT)-1) interior(iT).g_ms2 thisgcore];
    P_Planet_MPa(iT,:) = [P_MPa(iT,1:indSil(iT)-1) interior(iT).Pmantle_MPa thisPcore];
    T_Planet_K(iT,:) = [T_K(iT,1:indSil(iT)-1) interior(iT).Tmantle_K thisTcore];
    r_Planet_m(iT,:) = [r_m(iT,1:indSil(iT)-1)    interior(iT).r_mantle_m thisrcore];
    rho_pPlanet_kgm3(iT,:) = [rho_kgm3(iT,1:indSil(iT)-1) interior(iT).rho_mantle_kgm3 thisrhocore];
    ocean_thk(iT) = r_Planet_m(iT,indsLiquid(1)) - r_Planet_m(iT,indsLiquid(end)+1);
    ind_Obot(iT) = indsLiquid(end);
    
    VP_Planet_kms(iT,indsI) = velsIce.VIl_kms(iT,indsI);
    VP_Planet_kms(iT,indsLiquid) = vfluid_kms(iT,indsLiquid);
    VP_Planet_kms(iT,indsIII) = velsIce.VIIIl_kms(iT,indsIII);
    VP_Planet_kms(iT,indsV) = velsIce.VVl_kms(iT,indsV);
    VP_Planet_kms(iT,indsClath) = velsIce.Vclathl_kms(iT,indsClath);
    VP_Planet_kms(iT,indsVI) = velsIce.VVIl_kms(iT,indsVI);
    VP_Planet_kms(iT,indSil(iT):indSil(iT)+length(interior(iT).VP_mantle_kms)-1) = interior(iT).VP_mantle_kms;
    if Planet.FeCore
        start_core = indSil(iT)+length(interior(iT).VP_mantle_kms);
        VP_Planet_kms(iT,start_core:start_core-1+length(thisVPcore)) = thisVPcore;
    end
    
%     VS_Planet_kms(iT,:) = [velsIce.VIt_kms(iT,1:n_iceI) 0*vfluid_kms(iT,indsLiquid) ...
%         velsIce.VIIIt_kms(iT,indsIII) velsIce.VVt_kms(iT,indsV) velsIce.VVIt_kms(iT,indsVI) ...
%         interior(iT).VS_mantle_kms thisVScore];
    VS_Planet_kms(iT,indsI) = velsIce.VIt_kms(iT,indsI);
    VS_Planet_kms(iT,indsLiquid) = 0*vfluid_kms(iT,indsLiquid);
    VS_Planet_kms(iT,indsIII) = velsIce.VIIIt_kms(iT,indsIII);
    VS_Planet_kms(iT,indsV) = velsIce.VVt_kms(iT,indsV);
    VS_Planet_kms(iT,indsClath) = velsIce.Vclatht_kms(iT,indsClath);
    VS_Planet_kms(iT,indsVI) = velsIce.VVIt_kms(iT,indsVI);
    VS_Planet_kms(iT,indSil(iT):indSil(iT)+length(interior(iT).VS_mantle_kms)-1) = interior(iT).VS_mantle_kms;
    if Planet.FeCore
        VS_Planet_kms(iT,start_core:start_core-1+length(thisVPcore)) = thisVScore;
    end
    
    if isfield(mantle,'Ks_fn')
    %     Ks_Planet_GPa(iT,:)=[velsIce.KsI_GPa(iT,1:n_iceI) Ksfluid_GPa(iT,indsLiquid) ...
    %         velsIce.KsIII_GPa(iT,indsIII) velsIce.KsV_GPa(iT,indsV) velsIce.KsVI_GPa(iT,indsVI) ...
    %         interior(iT).Ks_GPa thisKScore];
        Ks_Planet_GPa(iT,indsI) = velsIce.KsI_GPa(iT,indsI);
        Ks_Planet_GPa(iT,indsLiquid) = Ksfluid_GPa(iT,indsLiquid);
        Ks_Planet_GPa(iT,indsIII) = velsIce.KsIII_GPa(iT,indsIII);
        Ks_Planet_GPa(iT,indsV) = velsIce.KsV_GPa(iT,indsV);
        Ks_Planet_GPa(iT,indsClath) = velsIce.Ksclath_GPa(iT,indsClath);
        Ks_Planet_GPa(iT,indsVI) = velsIce.KsVI_GPa(iT,indsVI);
        Ks_Planet_GPa(iT,indSil(iT):indSil(iT)+length(interior(iT).VS_mantle_kms)-1) = interior(iT).Ks_GPa;
        if Planet.FeCore
            Ks_Planet_GPa(iT,start_core:start_core-1+length(thisVPcore)) = thisKScore;
        end

    %     Gs_Planet_GPa(iT,:)=[velsIce.GsI_GPa(iT,1:n_iceI) 0*Ksfluid_GPa(iT,indsLiquid) ...
    %         velsIce.GsIII_GPa(iT,indsIII) velsIce.GsV_GPa(iT,indsV) velsIce.GsVI_GPa(iT,indsVI) ...
    %         interior(iT).Gs_GPa thisGScore];
        Gs_Planet_GPa(iT,indsI) = velsIce.GsI_GPa(iT,indsI);
        Gs_Planet_GPa(iT,indsLiquid) = 0*Ksfluid_GPa(iT,indsLiquid);
        Gs_Planet_GPa(iT,indsIII) = velsIce.GsIII_GPa(iT,indsIII);
        Gs_Planet_GPa(iT,indsV) = velsIce.GsV_GPa(iT,indsV);
        Gs_Planet_GPa(iT,indsClath) = velsIce.Gsclath_GPa(iT,indsClath);
        Gs_Planet_GPa(iT,indsVI) = velsIce.GsVI_GPa(iT,indsVI);
        Gs_Planet_GPa(iT,indSil(iT):indSil(iT)+length(interior(iT).VS_mantle_kms)-1) = interior(iT).Gs_GPa;
        if Planet.FeCore
            Gs_Planet_GPa(iT,start_core:start_core-1+length(thisVPcore)) = thisGScore;
        end

    end
%     QS_overfgamma_Planet(iT,:) = [QS_overfgamma_iceI(iT,:) 0*vfluid_kms(iT,indsLiquid)...
%         QS_overfgamma_iceIII(iT,indsIII) QS_overfgamma_iceV(iT,indsV) QS_overfgamma_iceVI(iT,indsVI) ...
%         interior(iT).QS_overfgamma thisQScore];

%% attenuation in ice
    Hp_iceI = Seismic.g_aniso_iceI*Planet.Tb_K(iT);
    QS_overfgamma_iceI = Seismic.B_aniso_iceI*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceI./T_K(iT,indsI))/Seismic.LOW_ICE_Q;
    try
        Tthis=T_K(iT,indsClath);
        Hp_clath = Seismic.g_aniso_clath*Planet.Tb_K(iT);
        QS_overfgamma_clath= Seismic.B_aniso_clath*....
            exp(Seismic.gamma_aniso_clath*Hp_clath./Tthis)/Seismic.LOW_ICE_Q;
    catch
        Hp_clath = Seismic.g_aniso_iceI*Planet.Tb_K(iT);
        QS_overfgamma_clath = Seismic.B_aniso_iceI*....
            exp(Seismic.gamma_aniso_iceI*Hp_clath./T_K(iT,indsClath))/Seismic.LOW_ICE_Q;
    end
    Tthis=T_K(iT,indsIII);
    Hp_iceIII = Seismic.g_aniso_iceIII*max(T_K(iT,indsIII));
    QS_overfgamma_iceIII = Seismic.B_aniso_iceIII*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceIII./Tthis)/Seismic.LOW_ICE_Q;
    Tthis=T_K(iT,indsV);
    Hp_iceV = Seismic.g_aniso_iceV*max(T_K(iT,indsV));
    QS_overfgamma_iceV = Seismic.B_aniso_iceV*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceV./Tthis)/Seismic.LOW_ICE_Q;
    Tthis=T_K(iT,indsVI);
    Hp_iceVI = Seismic.g_aniso_iceVI*max(T_K(iT,indsVI));
    QS_overfgamma_iceVI = Seismic.B_aniso_iceVI*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceVI./Tthis)/Seismic.LOW_ICE_Q;

    QS_overfgamma_Planet(iT,indsI) = QS_overfgamma_iceI;
    QS_overfgamma_Planet(iT,indsLiquid) = 0*vfluid_kms(iT,indsLiquid);
    QS_overfgamma_Planet(iT,indsIII) = QS_overfgamma_iceIII;
    QS_overfgamma_Planet(iT,indsV) = QS_overfgamma_iceV;
    QS_overfgamma_Planet(iT,indsVI) = QS_overfgamma_iceVI;
    QS_overfgamma_Planet(iT,indsClath) = QS_overfgamma_clath;
    QS_overfgamma_Planet(iT,indSil(iT):indSil(iT)+length(interior(iT).VS_mantle_kms)-1) = interior(iT).QS_overfgamma;
    if Planet.FeCore
        QS_overfgamma_Planet(iT,start_core:start_core-1+length(thisVPcore)) = thisQScore;
    end
    
 %%save the data to a text file
    Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [P_Planet_MPa(iT,:)', T_Planet_K(iT,:)', r_Planet_m(iT,:)'*1e-3, rho_pPlanet_kgm3(iT,:)',VP_Planet_kms(iT,:)',VS_Planet_kms(iT,:)',QS_overfgamma_Planet(iT,:)' Ks_Planet_GPa(iT,:)' Gs_Planet_GPa(iT,:)' g_Planet_ms2(iT,:)' phasePlanet(iT,:)'];
    header = sprintf('%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',...
        'P (MPa)','T (K)','r (km)','rho (kg m-3)','VP (km s-1)','VS (km s-1)','QS/gamma','KS (GPa)','GS (GPa)','g (m/s2)','phase');
    if cfg.CONDUCT %11/8/19 added g and phase. This changes the column number of electrical conductivity
        k_S_mPlanet(iT,indsI) = 0*QS_overfgamma_iceI;
        k_S_mPlanet(iT,indsLiquid) = k_S_m(iT,indsLiquid);
        k_S_mPlanet(iT,indsIII) = 0*QS_overfgamma_iceIII;
        k_S_mPlanet(iT,indsV) = 0*QS_overfgamma_iceV;
        k_S_mPlanet(iT,indsClath) = 0*QS_overfgamma_clath;
        k_S_mPlanet(iT,indsVI) = 0*QS_overfgamma_iceVI;
        k_S_mPlanet(iT,indSil(iT):indSil(iT)+length(interior(iT).VS_mantle_kms)-1) = 0*interior(iT).QS_overfgamma;
        if Planet.FeCore
            k_S_mPlanet(iT,start_core:start_core-1+length(thisVPcore)) = 0*thisQScore;
        end
        Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma k_S_mPlanet(iT,:)' ]; %#ok<AGROW>
        header = sprintf('%s%s\t%s\t',header,'k (S m-1)');
        
        kmean(iT) = mean(k_S_mPlanet(iT,indsLiquid));
        ktop(iT) = k_S_mPlanet(iT,indsLiquid(1));
    end
    if isfield(Planet,'Clathrate'); clathID = ['_Zclath' num2str(Zclath(iT)./1000,'%2.0f') 'km']; else; clathID = ''; end
    thissavestr = [savefile '_Zb' strLow num2str(Planet.zb_outerIce_m(iT)./1000,'%2.0f') 'km' clathID ];
    dlmwrite(fullfile([datpath thissavestr '_pp' num2str(iT) '.txt']),header,'delimiter','');
    dlmwrite(fullfile([datpath thissavestr '_pp' num2str(iT) '.txt']),Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma,...
        'delimiter','\t',...
        'precision','%3.5e',...
        '-append');
    dlmwrite(fullfile([datpath thissavestr '_mantle_pp' num2str(iT) '.dat']),[10*interior(iT).Pmantle_MPa' interior(iT).Tmantle_K'],'delimiter','\t');

    if isfield(Params,'TauP')
        taup_model=create_Taup([datpath thissavestr '.txt']);
        
    end
    
    if ~cfg.SKIP_PROFILES
    %% plot the seismic data and attenuation
        if Params.HOLD
            set(0, 'CurrentFigure', figs.seis);
            set(gcf,'Position', [476   662   857   505],'Name',lbl.seism)
        else
            set(0, 'CurrentFigure', figs.seis(iT));
            set(gcf,'Position', [476   662   857   505],'Name',[lbl.seism ' Tb = ' num2str(Planet.Tb_K(iT))])
            clf;
        end

        hp = subplot(1,3,1);
        if Params.HOLD
            hold on
        end
        plot(VS_Planet_kms(iT,:)',r_Planet_m(iT,:)'*1e-3,...
            VP_Planet_kms(iT,:)',r_Planet_m(iT,:)'*1e-3,'--','LineWidth',cfg.LW_seism)
            set(gcf,'color','w')
            xlabel('Sound Speeds (km s^{-1})','Fontsize',lbl.mltext)
            ylabel([math 'r_{' nm Planet.name '}' nm ' (km)'],'Fontsize',lbl.mltext);
            set(gca,'xlim',[0 1.1*max(VP_Planet_kms(iT,:))],'ylim',[0 ymax])
            grid on


        hp = subplot(1,3,2);
        if Params.HOLD
            hold on
        end
        if max(P_Planet_MPa(iT,:))<100
            plot(10*P_Planet_MPa(iT,:)',r_Planet_m(iT,:)'*1e-3,...
                T_Planet_K(iT,:)',r_Planet_m(iT,:)'*1e-3,'--',...
                rho_pPlanet_kgm3(iT,:)',r_Planet_m(iT,:)'*1e-3,'k-.',...
            'LineWidth',cfg.LW_seism);
            xlabel([math 'P' nm ' (bar), ' math 'T' nm ' (K), ' math '\rho' nm ' (kg m^{-3})'],'Fontsize',lbl.mltext);
            xmax = max([max(rho_pPlanet_kgm3(iT,:)) max(10*P_Planet_MPa(iT,:)) max(T_Planet_K(iT,:))]);
        else
            plot(P_Planet_MPa(iT,:)',r_Planet_m(iT,:)'*1e-3,...
                T_Planet_K(iT,:)',r_Planet_m(iT,:)'*1e-3,'--',...
                rho_pPlanet_kgm3(iT,:)',r_Planet_m(iT,:)'*1e-3,'k-.',...
            'LineWidth',cfg.LW_seism);
            xlabel([math 'P' nm ' (MPa), ' math 'T' nm ' (K), ' math '\rho' nm ' (kg m^{-3})'],'Fontsize',lbl.mltext);
            xmax = max([max(rho_pPlanet_kgm3(iT,:)) max(P_Planet_MPa(iT,:)) max(T_Planet_K(iT,:))]);
        end
        set(gca,'ylim',[0 ymax],'xlim',[0 xmax]);
        grid on

        subplot(1,3,3)
        if Params.HOLD
            hold on
        end
        hp = plot(QS_overfgamma_Planet(iT,:)',r_Planet_m(iT,:)'*1e-3,'LineWidth',cfg.LW_seism);
        set(gca,'xscale','log','ylim',[0 ymax],'xlim',[10 1e7],'XTick',[10 100 1000 1e4 1e5 1e6 1e7])
        grid on
        xlabel([math 'Q_S/\omega^{\gamma}'],'Fontsize',lbl.mltext)
        set(gcf,'color','w')

        if Params.HOLD
            set(0, 'CurrentFigure', figs.gsks);
            set(gcf,'Position', [476   662   857   505],'Name',lbl.gs_ks)
        else
            set(0, 'CurrentFigure', figs.gsks(iT));
            set(gcf,'Position', [476   662   857   505],'Name',[lbl.gs_ks ' Tb = ' num2str(Planet.Tb_K(iT))])
            clf;
        end
        hold on
        hp = plot(Ks_Planet_GPa(iT,:)',r_Planet_m(iT,:)'*1e-3,...
            Gs_Planet_GPa(iT,:)',r_Planet_m(iT,:)'*1e-3,'--','LineWidth',cfg.LW_seism);
        set(hp,'LineWidth',cfg.LW_seism);
        ylabel([math 'r_{' nm Planet.name '}' nm ' (km)'],'Fontsize',lbl.mltext);
        xlabel([math 'G_S' nm ' and ' math 'K_S' nm ' (GPa)'],'Fontsize',lbl.mltext);
        axis tight

        if ~cfg.HOLD
            
            print(figs.seis(iT),Params.savefigformat,fullfile([figpath thissavestr '_' vseis cfg.xtn]));
            print(figs.gsks(iT),Params.savefigformat,fullfile([figpath thissavestr '_' vgsks cfg.xtn]));
            if POROUS % Saved from earlier -- now we have thissavestr.
                print(figs.porP(iT),cfg.fig_fmt,fullfile([figpath thissavestr '_' vsP   cfg.xtn]));
                print(figs.porR(iT),cfg.fig_fmt,fullfile([figpath thissavestr '_' vsR   cfg.xtn]));
                print(figs.perm(iT),cfg.fig_fmt,fullfile([figpath thissavestr '_' vperm cfg.xtn]));
            end
        end
    end % ~cfg.SKIP_PROFILES
end

% Add conductivities and boundaries to Planet for usage in
% LayeredInductionResponse (but don't overwrite if we did something
% special with cfg.REDUCED)
if ~cfg.REDUCED
    Planet.sig = flip(k_S_mPlanet,2);
    Planet.boundaries = flip(r_Planet_m,2);
end

% We now have everything we need if we just want to do layered
% induction calculations.
% Everything left to do in PlanetProfile is for plotting purposes.
if ~cfg.SKIP_PROFILES
    
    % Save plots special if we are overlaying
    if cfg.HOLD
        print(figs.seis,Params.savefigformat,fullfile([figpath savebase vseis cfg.xtn]));
        print(figs.gsks,Params.savefigformat,fullfile([figpath savebase vgsks cfg.xtn]));
        if POROUS
            print(figs.porP,Params.savefigformat,fullfile([figpath savebase vsP   cfg.xtn]));
            print(figs.porR,Params.savefigformat,fullfile([figpath savebase vsR   cfg.xtn]));
            print(figs.perm,Params.savefigformat,fullfile([figpath savebase vperm cfg.xtn]));
        end
    end

    revcop = colormap('copper');
    revcop = revcop(256:-1:1,:);
    colormap(revcop);
    inferno_data = CM_inferno;
    colormap(inferno_data);

    set(0, 'CurrentFigure', figs.pvt6);
    if ~cfg.HOLD; clf; end
    clear opts
    set(gcf,'Position', [496   642   824   678],'Name',[lbl.inter ' x 6'])
    opts.Punits = 'MPa';
    opts.Ttight = true;
    if Planet.FeCore
        pinds = find(mantle.p(:,1)*1e3>=thisPcore(1));
        pcore = mantle.p(pinds,:)*1e3; tcore = mantle.t(pinds,:);
        if ~exist('core','var')
            core=deal([]);
        end
    else
        [core,tcore,pcore]=deal([]);
    end
    subplot(2,3,1);
    plotSolidInterior('rho','Density (kg m^{-3})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
    subplot(2,3,2); plotSolidInterior('cp',[math 'C_p' bnm ' (J m^{-3} K^{-1})'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
    subplot(2,3,3); plotSolidInterior('alpha',[math '\alpha' bnm ' (K^{-1})'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
    subplot(2,3,4); plotSolidInterior('vp',[math 'V_P' bnm ' (km s^{-1})'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
    subplot(2,3,5); plotSolidInterior('vs',[math 'V_S' bnm ' (km s^{-1})'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)

    if isfield(mantle,'Ks')
        subplot(2,3,6); plotSolidInterior('Ks',[math 'K_S' bnm ' (GPa)'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)

        set(0, 'CurrentFigure', figs.pvt4);
        if ~cfg.HOLD; clf; end
        set(gcf,'Position', [426   542   549   678],'Name',[lbl.inter ' x 4'])
        opts.Punits = 'GPa';
        colormap(inferno_data);
            subplot(2,2,1); plotSolidInterior('rho','Density (kg m^{-3})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
            subplot(2,2,2); plotSolidInterior('cp',[math 'C_p' bnm ' (J m^{-3} K^{-1})'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
            subplot(2,2,3); plotSolidInterior('Ks',[math 'K_S' bnm ' (GPa)'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
            subplot(2,2,4); plotSolidInterior('Gs',[math 'G_S' bnm ' (GPa)'],T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)

        print(figs.pvt4,Params.savefigformat,fullfile([figpath savefile '_' vpvt4 cfg.xtn]));
    end

    % save the mineral compositions plot
    print(figs.pvt6,Params.savefigformat,fullfile([figpath savefile '_' vpvt6 cfg.xtn]));

    %% create the legend that describes tb, z_b, and heat flux
    lstr_2 = {};
    for iT = 1:nTbs
        lstr_2 = {lstr_2{:} [math 'T_{b}' nm ': ' num2str(Planet.Tb_K(iT),'%0.2f') ' K, ' math 'q_{b}/q_{c}' nm ':'...
            num2str(1e3*Qb(iT),'%0.0f') '/' num2str(1e3*Q_Wm2(iT),'%0.0f') ' mW m^{-2}, ' math 'z_{b}' nm ':' num2str(1e-3*Planet.zb_outerIce_m(iT),'%0.0f') ' km']};
    end

    set(0, 'CurrentFigure', figs.wedg);
    clf;
    for iT = 1:nTbs
        subplot(1,nTbs,iT); hold on
        input.phase = phasePlanet(iT,:);
        input.r_m = r_Planet_m(iT,:);
        PlanetWedge(input);
        box off
        if iT==1
            ylabel(Planet.name,'FontSize',24)
        else
          set(gcf, 'Name', lbl.wedge)
          set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none')

        end
        title([math 'T_b' nm ' = ' num2str(Planet.Tb_K(iT)) ' K'],'FontSize',lbl.lgtext)
    end
    print(figs.wedg,Params.savefigformat,fullfile([figpath savefile '_' vwedg cfg.xtn]));


    %%  plot profile with 4 subplots
    %%
    if Params.foursubplots
    Dsil_km=(Planet.R_m-R_sil_mean_m(:))*1e-3;

    set(0, 'CurrentFigure', figs.cond);
    maxScale = 1.01;
    minScale = 0.99;
    if ~Params.HOLD
        clf;
    end
        set(gcf,'Position', [285 33 898 751], 'Name', lbl.panl4)
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        subplot(3,6,4:6);
    else
        subplot(2,6,4:6);
    end
        hold on
    for iT = 1:nTbs
        line(T_K(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),...
            'LineWidth',cfg.LW_seism,'LineStyle',LineStyle);
        hm = line(T_K(iT,indSil(iT)),z_m(iT,indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'Marker','o');
        if strcmp(LineStyle,'-')
            hm.MarkerFaceColor=Params.colororder(:,iT);
        end
    end
    if Params.HOLD
        ax = gca;
        if max(max(T_K))>ax.XLim
            ax.XLim(2)=max(max(T_K));
        end
        if maxScale*max(Dsil_km)>ax.YLim
            ax.YLim(2) = maxScale*max(Dsil_km);
        end
    end
        set(gca,'ydir','reverse',...
            'xlim',[245 max(max(T_K))],'ylim',[0 maxScale*max(Dsil_km)],...
            'FontSize',lbl.mdtext,'YAxisLocation','right');%'XAxisLocation','top');
        set(gcf, 'Name', lbl.panl4);
    
    xlabel('Temperature (K)','Fontsize',lbl.mltext);
    ylabel('Depth (km)','Fontsize',lbl.mltext);
    % for iT=1:nTbs
    %     hline(Dsil_km(iT),Params.colororder(:,iT));
    % end
    box on;

    %== Plot sounds speeds in the ice and liquid layers
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        subplot(3,6,10);
    else
        subplot(2,6,10);
    end
    hold on;
    % clear hp; % putting this here because an error was raised when running
    % Triton models. May 4 2020. I think the problem is the model I'm running
    % doesn't have any liquid. yep.
    for iT = 1:nTbs
        hp(iT) =  line(vfluid_kms(iT,Params.nsteps_iceI+1:indSil(iT)),z_m(iT,Params.nsteps_iceI+1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        maxV(iT) = max(vfluid_kms(iT,Params.nsteps_iceI+1:indSil(iT)));
        minV(iT) = min(vfluid_kms(iT,Params.nsteps_iceI+1:indSil(iT)));
    end
    maxV = max(maxV); minV = min(minV);
        ax = gca;

    if Params.HOLD
        if max(Dsil_km)>ax.YLim
            ax.YLim(2) = max(Dsil_km);
        end
        if max(maxV)>ax.XLim(2)
            ax.XLim(2) = maxV;
        end
            if maxV>ax.XLim(2)
            ax.XLim(2) = maxV;
        end
        if minV<ax.XLim(1)
            ax.XLim(1) = minV;
        end
    else
       ax.XLim = [minV maxV];
       ax.YLim = [0 max(Dsil_km)];
    end
       ax.YDir = 'reverse';
       ax.FontSize = lbl.mdtext;

       box on

       if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        subplot(3,6,11);
        else
        subplot(2,6,11);
       end

    hold on;
    for iT = 1:nTbs
        line(velsIce.VIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VVt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VVIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
           line(velsIce.Vclatht_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
       maxV(iT) = max([velsIce.VIt_kms(iT,1:indSil(iT)) ...
            velsIce.VIIt_kms(iT,1:indSil(iT)) ...
            velsIce.Vclatht_kms(iT,1:indSil(iT))...
            velsIce.VIIIt_kms(iT,1:indSil(iT)) ...
            velsIce.VVt_kms(iT,1:indSil(iT)) ...
            velsIce.VVIt_kms(iT,1:indSil(iT))]);
       minV(iT) = min([velsIce.VIt_kms(iT,1:indSil(iT)) ...
         velsIce.Vclatht_kms(iT,1:indSil(iT)) ...
            velsIce.VIIt_kms(iT,1:indSil(iT)) ...
            velsIce.VIIIt_kms(iT,1:indSil(iT)) ...
            velsIce.VVt_kms(iT,1:indSil(iT)) ...
            velsIce.VVIt_kms(iT,1:indSil(iT))]);
    end
    maxV = max(maxV);minV = min(minV);

    ax = gca;
    if Params.HOLD
        if max(Dsil_km)>ax.YLim
            ax.YLim(2) = max(Dsil_km);
        end
        if maxV>ax.XLim(2)
            ax.XLim(2) = maxV;
        end
        if minV<ax.XLim(1)
            ax.XLim(1) = minV;
        end
    else
       ax.XLim = [minV maxV];
       ax.YLim = [0 max(Dsil_km)];
    end
       ax.YDir = 'reverse';
       ax.FontSize = lbl.mdtext;
       ax.YTickLabel = [];
    box on
    xlabel('Sound Speed (km s^{-1})','Fontsize',lbl.mltext);


    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        subplot(3,6,12);
    else
        subplot(2,6,12);
    end
    hold on;
    for iT = 1:nTbs
        line(velsIce.VIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VVl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
         line(velsIce.Vclathl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_sound);
    
    maxV(iT) = max([velsIce.VIl_kms(iT,1:indSil(iT)) ...
            velsIce.Vclathl_kms(iT,1:indSil(iT)) ...
            velsIce.VIIl_kms(iT,1:indSil(iT)) ...
            velsIce.VIIIl_kms(iT,1:indSil(iT)) ...
            velsIce.VVl_kms(iT,1:indSil(iT)) ...
            velsIce.VVIl_kms(iT,1:indSil(iT)) ...
            velsIce.Vclathl_kms(iT,1:indSil(iT))]);
    minV(iT) = min([velsIce.VIl_kms(iT,1:indSil(iT))  ...
            velsIce.Vclathl_kms(iT,1:indSil(iT)) ...
            velsIce.VIIl_kms(iT,1:indSil(iT)) ...
            velsIce.VIIIl_kms(iT,1:indSil(iT)) ...
            velsIce.VVl_kms(iT,1:indSil(iT)) ...
            velsIce.VVIl_kms(iT,1:indSil(iT))]);
    end
    maxV = max(maxV);minV = min(minV);

    ax = gca;
    if Params.HOLD
        if max(Dsil_km)>ax.YLim
            ax.YLim(2) = max(Dsil_km);
        end

        if maxV>ax.XLim(2)
            ax.XLim(2) = maxV;
        end
        if minV<ax.XLim(1)
            ax.XLim(1) = minV;
        end
    else
       ax.XLim = [minV maxV];
       ax.YLim = [0 max(Dsil_km)];
    end
       ax.YDir = 'reverse';
       ax.FontSize = lbl.mdtext;
          ax.YAxisLocation = 'right';


       ylabel('Depth (km)','Fontsize',lbl.mltext);

    box on

    % not sure this is used anymore. flagged for deletion
    % for iT=1:nTbs
    %     if ~isnan(velsIce.VVIl_kms(iT,indSil(iT)))
    %         velT = velsIce.VVIt_kms(iT,indSil(iT));
    %         velL = velsIce.VVIl_kms(iT,indSil(iT));
    %     elseif  ~isnan(velsIce.VVl_kms(iT,indSil(iT)))
    %         velT = velsIce.VVt_kms(iT,indSil(iT));
    %         velL = velsIce.VVl_kms(iT,indSil(iT));
    %     elseif  ~isnan(velsIce.VIIIl_kms(iT,indSil(iT)))
    %         velT = velsIce.VIIIt_kms(iT,indSil(iT));
    %         velL = velsIce.VIIIl_kms(iT,indSil(iT));
    %     elseif  ~isnan(velsIce.VIIl_kms(iT,indSil(iT)))
    %         velT = velsIce.VIIt_kms(iT,indSil(iT));
    %         velL = velsIce.VIIl_kms(iT,indSil(iT));
    %     elseif  ~isnan(vfluid_kms(iT,indSil(iT)))
    %         velT = vfluid_kms(iT,indSil(iT));
    %         velL = vfluid_kms(iT,indSil(iT));
    %     else
    %         velT = velsIce.VIt_kms(iT,indSil(iT));
    %         velL = velsIce.VIt_kms(iT,indSil(iT));
    %     end
    % end

    % ax1 = gca;
    % ax1.XAxisLocation = 'top';
    % ax1.YAxisLocation = 'right';
    % ax1_pos = ax1.Position;
    % ax2 = axes('Position',ax1_pos,...
    %     'XAxisLocation','bottom',...
    %     'YAxisLocation','left',...
    %     'YColor','none',...
    %     'Color','none');


     if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
       if wo>0
        subplot(3,6,16:18);hold on
        switch Planet.Ocean.comp
            case 'MgSO4'
                for iT = 1:nTbs
                    line(k_S_m(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,...
                        'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_std);
                    mink_val(iT) = min(k_S_m(iT,1:indSil(iT)));
                    maxk_val(iT) = max(k_S_m(iT,1:indSil(iT)));
                end
            case 'Seawater'
                for iT = 1:nTbs
                    line(k_S_m(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),...
                        'LineStyle',LineStyle,'LineWidth',cfg.LW_std);
                    mink_val(iT) = min(k_S_m(iT,1:indSil(iT)));
                    maxk_val(iT) = max(k_S_m(iT,1:indSil(iT)));
                end
            case 'NH3'
                for iT = 1:nTbs
    %                     line(k_S_m(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle',LineStyle,'LineWidth',cfg.LW_std);
                    mink_val(iT) = NaN;
                    maxk_val(iT) = NaN;
                end
            case 'NaCl'
                for iT = 1:nTbs
                    mink_val(iT) = NaN;
                    maxk_val(iT) = NaN;
                end
        end
        mink_val = min(mink_val);maxk_val = max(maxk_val);

        ax = gca;
        if Params.HOLD
            if max(Dsil_km)>ax.YLim(2)
                ax.YLim(2) = max(Dsil_km);
            end
            if maxk_val>ax.XLim(2)
                ax.XLim(2) = maxk_val;
            end
            if mink_val<ax.XLim(1)
                ax.XLim(1) = mink_val;
            end
        else
            if ~isnan(mink_val)
                ax.XLim = [mink_val maxk_val];
                ax.YLim = [0 max(Dsil_km)];
            end
        end
        ax.YDir = 'reverse';
        ax.FontSize = lbl.mdtext;
        ax.YAxisLocation = 'right';

        xlabel('Electrical Conductivity (S m^{-1})','Fontsize',lbl.mltext);
        ylabel('Depth (km)','Fontsize',lbl.mltext);
        box on
       end
    end


    % subplot(3,2,[1 3 5]);
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        subplot(3,6,[1:3 7:9 13:15]);
    else
        subplot(2,6,[1:3 7:9]);
    end
    hold on;
    % hold all
    R2ind = C2mean+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
    Psil_MPa = diag(P_MPa(:,R2ind(:)));
    for iT=1:nTbs
        ht(iT)= plot(P_MPa(iT,1:indSil(iT)),rho_kgm3(iT,1:indSil(iT)),'Color',Params.colororder(:,iT),...
            'LineWidth',cfg.LW_std,'LineStyle',LineStyle);
        hm = plot(Psil_MPa(iT),interp1(P_MPa(iT,1:end-2),rho_kgm3(iT,1:end-2),Psil_MPa(iT)),'Color',Params.colororder(:,iT));
        if ~Params.HOLD
            hm.MarkerFaceColor=Params.colororder(:,iT);
        end
    end

    for ir = 1:length(Params.wref)
        hw(ir) = plot(Pref_MPa,rho_ref_kgm3(ir,:),['k' Params.wrefLine]);
    end

    if Params.Legend
        hleg1 = legend([ht],lstr_2{:},'Fontsize',lbl.smtext);
        set(hleg1,'location',Params.LegendPosition,'box','off')
    end

    box on;
    ax = gca;
    ax.FontSize = lbl.mdtext;
     if Params.HOLD
        if maxScale*max(max((rho_kgm3(:,1:indSil(:)))))>ax.YLim(2)
            ax.YLim(2) = maxScale*max(max((rho_kgm3(:,1:indSil(:)))));
        end
        if maxScale*max(Psil_MPa)>ax.XLim(2)
            ax.XLim(2) = maxScale*max(Psil_MPa);
        end
     else
        ax.XLim = [0 maxScale*max(Psil_MPa)];
        ax.YLim = Params.ylim;
     end
    xlabel('Pressure (MPa)','Fontsize',lbl.mltext)
    ylabel('Density (kg m^{-3})','Fontsize',lbl.mltext)

    print(figs.cond,Params.savefigformat,fullfile([figpath savebase vcond cfg.xtn]));

    %% plot profile with 2 subplots
    else

    set(0, 'CurrentFigure', figs.cond);
    maxScale = 1.01;
    if ~Params.HOLD
        clf;
    end
        set(gcf,'Position', [411    52   898    751], 'Name', lbl.panl4)
        subplot(1,2,2);
        hold on

    Dsil_km=(Planet.R_m-R_sil_mean_m(:))*1e-3;

    for iT = 1:nTbs
        hp(iT) =  line(vfluid_kms(iT,Params.nsteps_iceI+1:indSil(iT)),z_m(iT,Params.nsteps_iceI+1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',cfg.LW_sound);
     line(velsIce.Vclathl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.Vclatht_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','-.','LineWidth',SoundLineWidth);

    end
    for iT = 1:nTbs
        if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
            line(k_S_mMgSO4p01Planet(iT,1:indSil(iT))*25,z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineStyle','--','LineWidth',cfg.LW_seism);
        end
    end
    set(gca,'ydir','reverse','xlim',[0 5],'ylim',[0 maxScale*max(Dsil_km)]);%,'xlim',[1 4]);
        if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
            xlabel('Conductivity (S m^{-1} \times 25) and Sound Speed (km s^{-1})','Fontsize',lbl.mltext);
        else
            xlabel('Sound Speed (km s^{-1})','Fontsize',lbl.mltext);
        end
        ylabel('Depth (km)');
    box on

    % not sure this is used anymore. flagged for deletion
    %
    % for iT=1:nTbs
    %     if ~isnan(velsIce.VVIl_kms(iT,indSil(iT)))
    %         velT = velsIce.VVIt_kms(iT,indSil(iT));
    %         velL = velsIce.VVIl_kms(iT,indSil(iT));
    %     elseif  ~isnan(velsIce.VVl_kms(iT,indSil(iT)))
    %         velT = velsIce.VVt_kms(iT,indSil(iT));
    %         velL = velsIce.VVl_kms(iT,indSil(iT));
    %     elseif  ~isnan(velsIce.VIIIl_kms(iT,indSil(iT)))
    %         velT = velsIce.VIIIt_kms(iT,indSil(iT));
    %         velL = velsIce.VIIIl_kms(iT,indSil(iT));
    %     elseif  ~isnan(velsIce.VIIl_kms(iT,indSil(iT)))
    %         velT = velsIce.VIIt_kms(iT,indSil(iT));
    %         velL = velsIce.VIIl_kms(iT,indSil(iT));
    %     elseif  ~isnan(vfluid_kms(iT,indSil(iT)))
    %         velT = vfluid_kms(iT,indSil(iT));
    %         velL = vfluid_kms(iT,indSil(iT));
    %     else
    %         velT = velsIce.VIt_kms(iT,indSil(iT));
    %         velL = velsIce.VIt_kms(iT,indSil(iT));
    %     end
    % end

    ax1 = gca;
    ax1.XAxisLocation = 'top';
    ax1.YAxisLocation = 'right';
    ax1_pos = ax1.Position;
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','bottom',...
        'YAxisLocation','left',...
        'YColor','none',...
        'Color','none');

    for iT = 1:nTbs
        line(ax2,T_K(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'LineWidth',cfg.LW_seism);
        hm = line(T_K(iT,indSil(iT)),z_m(iT,indSil(iT))*1e-3,'Color',Params.colororder(:,iT),'Marker','o');
        hm.MarkerFaceColor=Params.colororder(:,iT);
    end

    set(gca,'ydir','reverse','xlim',[245 320],'ylim',[0 maxScale*max(Dsil_km)]);
    xlabel('Temperature (K)','Fontsize',lbl.mltext);
    ylabel('Depth (km)','Fontsize',lbl.mltext);
    % for iT=1:nTbs
    %     hline(Dsil_km(iT),Params.colororder(:,iT));
    % end
    box on;

    subplot(1,2,1);
    hold all
    R2ind = C2mean+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
    Psil_MPa = diag(P_MPa(:,R2ind(:)));
    for iT=1:nTbs
        ht(iT)=  plot(P_MPa(iT,1:indSil(iT)),rho_kgm3(iT,1:indSil(iT)),'Color',Params.colororder(:,iT),'LineWidth',cfg.LW_seism);
        hm = plot(Psil_MPa(iT),interp1(P_MPa(iT,:),rho_kgm3(iT,:),Psil_MPa(iT)),'Color',Params.colororder(:,iT));
        hm.MarkerFaceColor=Params.colororder(:,iT);
    end

    hw(1) = plot(Pref_MPa,rho_ref_kgm3(1,:),'k--');
    hw(2) = plot(Pref_MPa,rho_ref_kgm3(2,:),'k--');
    hw(3) = plot(Pref_MPa,rho_ref_kgm3(3,:),'k--');
    hw(4) = plot(Pref_MPa,rho_ref_kgm3(4,:),'k--');


    % hleg1 = legend([ht hw],{lstr_2{:},'0 wt%','5 wt%','10 wt%','15 wt%'});%,'20 wt%');
    if Params.Legend
        hleg1 = legend([ht],lstr_2{:},'Fontsize',lbl.smtext);
        set(hleg1,'location','southeast','box','off')
    end

    axis tight; box on;
    set(gca,'xlim',[0 maxScale*max(Psil_MPa)])
    xlabel('Pressure (MPa)','Fontsize',lbl.mltext)
    ylabel('Density (kg m^{-3})','Fontsize',lbl.mltext)

    % ax(1) = gca;
    % ax(2)=axes('Position',get(ax(1),'Position'),...
    %    'XAxisLocation','top',...
    %    'YAxisLocation','right',...
    %    'XColor','none',...
    %    'Color','none');
    %
    % % add the sound speeds
    % for iT = 1:nTbs
    %     line(P_MPa(iT,:),vfluid_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VIl_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VIt_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VIIl_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VIIt_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VIIIl_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VIIIt_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VVl_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VVt_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VVIl_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    %     line(P_MPa(iT,:),velsIce.VVIt_kms(iT,:),'Color',Params.colororder(:,iT),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
    % end
    % set(ax(2),'ylim',[1.5 5])
    % set(ax(2),'xlim',[0 1800])
    % ylabel('Sound Speeds (km s^{-1})')

        print(figs.cond,Params.savefigformat,fullfile([figpath savefile cfg.xtn]));
    end
    %%

    set(0, 'CurrentFigure', figs.grav);
    set(gcf, 'Name', lbl.gravg)
    if ~Params.HOLD
        clf;
    end
    for iT = 1:nTbs
        gp(iT) = g_ms2(iT,C2mean(iT));
    end

    subplot(1,2,1);
    if Params.HOLD
        hold on;
    end
    % plot(g_ms2',r_m'*1e-3,gp,R_sil_mean_m*1e-3,'o')
    hl = plot(g_Planet_ms2',r_Planet_m'*1e-3);
    hold on
    for iT = 1:nTbs
        hp = plot(gp(iT),R_sil_mean_m*1e-3,'o');
        set(hl(iT),'LineWidth',cfg.LW_seism,'LineStyle',LineStyle,'Color',Params.colororder(:,iT));
        set(hp,'Color',Params.colororder(:,iT));
    end
    xlabel('g (m s^{-2})','Fontsize',lbl.mltext)
    ylabel('R (km)','Fontsize',lbl.mltext)
    axis tight

    subplot(1,2,2);
    if Params.HOLD
        hold on;
    end
    hl = plot(P_Planet_MPa',r_Planet_m'*1e-3);
    for iT = 1:nTbs
        set(hl(iT),'LineWidth',cfg.LW_seism,'LineStyle',LineStyle,'Color',Params.colororder(:,iT));
    end

     xlabel('Pressure (MPa)','Fontsize',lbl.mltext)
     ylabel([math 'r_{' nm Planet.name '}' nm ' (km)'],'Fontsize',lbl.mltext)
    axis tight

    print(figs.grav,Params.savefigformat,fullfile([figpath savebase vgrav cfg.xtn]));

    %%
    % figure(232);clf;
    % subplot(2,1,2)
    % hold all;
    % clear ht hw
    % for iT=1:nTbs
    %     ht(iT)=  plot(P_MPa(iT,:),T_K(iT,:),Params.colororder(:,iT));
    % end
    % for iw = 1:length(Tref_K(:,1))
    %     hw(iw) = plot(Pref_MPa,Tref_K(iw,:),'k--');
    % end
    %
    %
    % hleg2 = legend([ht hw],{lstr_2{:},'0 wt%','5 wt%','10 wt%','15 wt%'});
    % set(hleg2,'location','NorthEast')
    %
    % axis tight; box on;
    % %     title('Temperature as a function of Pressure in Planet for T_{s} = 110 K, T_{b} = 260K')
    % xlabel('Pressure (MPa)')
    % ylabel('Temperature (K)')
    fprintf(['\n @@@  Call complete, figures printed to ' figpath ' @@@\n\n'])
end % ~cfg.SKIP_PROFILES

outPlanet = Planet;
end %PlanetProfile
%% properties of water ice
%Note: for VspChoukroun inds: liquid = 1, Ice I = 2, II = 3, III = 4, V = 5, VI = 6
function rho_kgm3 = getRhoIce(P_MPa,T_K,ind)

for jj=1:length(ind)
    switch ind(jj)
        % 1 (liquid water) isn't possible because of the above convention used for
        % implementing Choukroun and Grasset's (2010) EOS.
        % MJS 2020-10-16: I don't think the above comment is applicable.
        % We can adjust the indices as needed, this is not a model limitation.
        case 0
            material = 'water';
            warning('WARNING: getRhoIce was called for liquid water.')
        case 1
            material = 'Ih';
        case 2
            material = 'II';
        case 3
            material = 'III';
        case 5
            material = 'V';
        case 6
            material = 'VI';
        case 30
            material = 'Clath';
            % clath_out = Helgerud_sI(P_MPa(jj),T_K(jj));
            %rho_kgm3(jj)=clath_out.rho;
    end
    try
        out = SeaFreeze([P_MPa(jj),T_K(jj)],material);
        rho_kgm3(jj) = out.rho;
    catch
        if ind(jj)<3 % ices I and II
            rho_kgm3(jj) = 1000./getVspChoukroun2010(P_MPa(jj),T_K(jj),ind(jj)+1);
        else
            MW_H2O = 18.01528;
            switch ind(jj)
                case 3 % ice III
                    iceIII_Vfn = load('iceIII_sp_V_fPT');
                    rho_kgm3(jj) = 1./fnval(iceIII_Vfn.sp_V_fPT,{1e-3*P_MPa(jj) T_K(jj)});
                case 5 % ice V
                    iceV_Vfn = load('iceV_sp_V_fPT');
                    rho_kgm3(jj) = 1./fnval(iceV_Vfn.sp_V_fPT,{1e-3*P_MPa(jj) T_K(jj)});
                case 6 % ice VI
                    rho_kgm3(jj) = MW_H2O*1000./iceVI_PT_EOS_bezacier(1e-3*P_MPa(jj),T_K(jj));
                case 30
                    clath_out = Helgerud_sI(P_MPa(jj),T_K(jj));
                    rho_kgm3(jj)=clath_out.rho;
            end
        end
    end
end

end %getRhoIce
function [Cp, alpha]= getCpIce(P_MPa,T_K,ind)
    % convert phase from our nomenclature to array index for
    % CpH2O_Choukroun
    if ind<=3
        ind = ind+1;
    end
    switch ind
        % 1 (liquid water) isn't possible because of the above convention used for
        % implementing Choukroun and Grasset's (2010) EOS.
        % MJS 2020-10-16: I don't think the above comment is applicable.
        % We can adjust the indices as needed, this is not a model limitation.
        case 1
            material = 'water';
            warning('WARNING: getRhoIce was called for liquid water.')
        case 2
            material = 'Ih';
        case 3
            material = 'II';
        case 4
            material = 'III';
        case 5
            material = 'V';
        case 6
            material = 'VI';
        case 30
            % convert to array index for use in CpH2O_Choukroun
            ind = 7;
            material = 'Clath';
    end
    try
        out = SeaFreeze([P_MPa,T_K],material);
        Cp = out.Cp;
        alpha = out.alpha;
    catch
        disp(['T_ice = ' num2str(T_K) '. This seems to be too low for SeaFreeze. Using Choukroun and Grasset (2010) instead.']);
        Cp = CpH2O_Choukroun(P_MPa,T_K,ind);
         cpout=Helgerud_sI(P_MPa,T_K); % Danger! Cp for ice
 
          alpha=cpout.alpha; % use better alpha
    end
end %getCpIce
function phase = getIcePhase(P_MPa,T_K,w_pct,str_comp)
    switch str_comp
        case 'MgSO4'
%             if P_MPa<800
%                 W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
%                 mo=1000./W_MgSO4./(1./(0.01.*w_pct)-1); %conversion to molality from wt%
%                 [~,phase] = TmeltMgSO4Bollengier(P_MPa,mo); % this is faster and more accurate
%             else
                phase = getIcePhaseMgSO4(P_MPa,T_K,w_pct);
                % these conditions only apply if the Bollengier function is
                % ignored
%                 if phase==2
%                     phase=1;
%                 elseif phase==4
%                     phase=3;
%                 end
%             end
        case 'NaCl'
%             phase = getIcePhaseMgSO4(P_MPa,T_K,w_pct);
            phase = LBFIcePhase(P_MPa,T_K,w_pct,'NaCl');
        case 'Seawater'
            global swEOS
            phase=swEOS.gsw.tfreezing(w_pct,10*P_MPa)>T_K;
        case 'NH3'
%             phase = getIcePhaseNH3(P_MPa,T_K,w_pct);
            phase = LBFIcePhase(P_MPa,T_K,w_pct,'NH3');
    end
end % getIcePhase
function Tfreeze_K = getTfreeze(P_MPa,wo,str_comp,Tprior)
    switch str_comp
        case 'MgSO4'
%             if P_MPa<800
%                 W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
%                 mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from wt%
%                 Tfreeze_K = TmeltMgSO4Bollengier(P_MPa,mo); % this is faster and more accurate
%             else
                Tfreeze_K = fzero(@(T) L_IceMgSO4(P_MPa,T,wo,1),[200 400]);
%             end
        case 'NaCl'
%             Pfreeze_MPa = fzero(@(P) 0.5-LBFIcePhase(P,T_K,wo,'NaCl'),[0.1 200]);
            options = optimset('fzero');
            options.TolX = 1e-7; % 1e-13 failed for 5wt% NaCl for Titan with Tb=258.0842 at P = 1163, step il = 280
            if isempty(Tprior)
                Trange = [240 400];
            else
                Trange = Tprior + [-1 1]*4;
            end
                %Tfreeze_K = fzero(@(T) 0.5-LBFIcePhase(P_MPa,T,wo,'NaCl'),Trange,options);
        case 'Seawater'
            global swEOS
            Tfreeze_K = swEOS.gsw.tfreezing(wo,10*P_MPa);
        case 'NH3'
%             Tfreeze_K = fzero(@(T) L_IceNH3(P_MPa,T,wo,1),[200 350]);
            options = optimset('fzero');
            options.TolX = 1e-7; % 1e-13 failed for 5wt% NaCl for Titan with Tb=258.0842 at P = 1163, step il = 280
            if isempty(Tprior)
                Trange = [239 401];
            else
                Trange = Tprior + [-1 1]*4;
            end
            try
                Tfreeze_K = fzero(@(T) 0.5-LBFIcePhase(P_MPa,T,wo,'NH3'),Trange,options);
            catch
                x = 1 ;
            end
    end
end % getTfreeze
function Pfreeze_MPa = getPfreeze(T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'
            Pfreeze_MPa = fzero(@(P) L_IceMgSO4(P,T_K,wo,1),[0 250]);
        case 'NaCl'
            options = optimset('fzero');
            options.TolX = 1e-11; % prevent errors from overthinking by fzero. the error is:
%             Operands to the || and && operators must be convertible to logical
% scalar values.
%
% Error in fzero (line 444)
% while fb ~= 0 && a ~= b
% it occurred while trying to find Pfreeze for 5wt% NaCl on Titan with Tb
% =258K and TolX = 1e-15
            %Pfreeze_MPa = fzero(@(P) 0.5-LBFIcePhase(P,T_K,wo,'NaCl'),[0.1 209],options); % the P bounds are sensitive because multivalued solutions (ices Ih and III) will cause fzero to fail
        case 'Seawater'
            global swEOS
            Pfreeze_MPa = 0.1*swEOS.gsw.pfreezing(wo,T_K);
        case 'NH3'
%             Pfreeze_MPa = fzero(@(P) L_IceNH3(P,T_K,wo,1),[0 250]);
            %LBF
            options = optimset('fzero');
            options.TolX = 1e-11; % prevent errors from overthinking by fzero. the error is:
            Pfreeze_MPa = fzero(@(P) 0.5-LBFIcePhase(P,T_K,wo,'NH3'),[0.1 209],options); % the P bounds are sensitive because multivalued solutions (ices Ih and III) will cause fzero to fail
    end
end %getPfreeze
function Pfreeze_MPa = getPfreezeIII(T_K,wo,str_comp)
    switch str_comp
       case 'MgSO4'
           Pfreeze_MPa = fzero(@(P) L_IceMgSO4(P,T_K,wo,3),[0 500]);
       case 'NaCl'
%            Pfreeze_MPa = fzero(@(P) 0.5-LBFIcePhase(P,T_K,wo,'NaCl'),[0.1 500]);
       case 'NH3'
           %Pfreeze_MPa = fzero(@(P) L_IceNH3(P,T_K,wo,3),[0 500]);
    end
end %getPfreezeIII
%% fluid properties
function [rho_kgm3,Cp,alpha_Km1]=fluidEOS(P_MPa,T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'
            W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
            mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from wt%
            [rho_kgm3,Cp,alpha_Km1]=MgSO4_EOS2_planetary_smaller(mo,P_MPa,T_K-273.15);
        case 'NaCl'
            error('ERROR: NaCl is not yet implemented for fluidEOS.')
        case 'Seawater'
            global swEOS
            extrapflag = 1;
            % expect that wo is absolute salinity
            rho_kgm3=swEOS.gsw.dens(wo,T_K,P_MPa*10);
            Cp = swEOS.gsw.cp(wo,T_K,P_MPa*10);
            alpha_Km1 = swEOS.gsw.alpha(wo,T_K,P_MPa*10);
            if isnan(rho_kgm3) && extrapflag
                    disp('WARNING: extrapolating fluid rho, Cp, and alpha; this is seawater above 120 MPa, right?')
                Pin = [110 115 120];
                rhoin=swEOS.gsw.dens(ones(1,3)*wo,ones(1,3)*T_K,Pin*10);
                Cpin = swEOS.gsw.cp(ones(1,3)*wo,ones(1,3)*T_K,Pin*10);
                alphain = swEOS.gsw.alpha(ones(1,3)*wo,ones(1,3)*T_K,Pin*10);
                
                rho_kgm3 = interp1(Pin,rhoin,P_MPa,'linear','extrap');
                Cp = interp1(Pin,Cpin,P_MPa,'linear','extrap');
                alpha_Km1 = interp1(Pin,alphain,P_MPa,'linear','extrap');
            end
        case 'NH3'
           %[rho_kgm3,~,~,Cp,~,alpha_Km1] = ...
           %     refproppy32([P_MPa T_K],{'ammonia' 'water'},[wo 100-wo]/100,-1);
%             error('ERROR: NH3 is not yet implemented for fluidEOS.')
            W_NH3 = 17.031;            
            mo=1000./W_NH3./(1./(0.01.*wo)-1); %conversion to molality from wt%
            out = SeaFreeze({P_MPa T_K mo},'NH3');
            rho_kgm3 = out.rho;
            Cp = out.Cp;
            alpha_Km1 = out.alpha;

    end
end %fluidEOS
function [vel_kms,Ks_GPa] = fluidSoundSpeeds(P_MPa,T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'
            W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
            mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from wt%
            vel_kms = MgSO4_EOS2_planetary_velocity_smaller_vector(mo,P_MPa,T_K-273.15); % assumes P and T are the same length.
        case 'NaCl'
            global swEOS
            disp('WARNING: NaCl is not yet implemented for fluidSoundSpeeds.')
        case 'Seawater'
            global swEOS
            % expect that wo is absolute salinity
            vel_kms=1e-3*swEOS.gsw.vel(wo*ones(1,length(T_K)),T_K,P_MPa*10);
        case 'NH3'
%             disp('WARNING: NH3 is not yet implemented for fluidSoundSpeeds.')
            W_NH3 = 17.031;            
            mo=1000./W_NH3./(1./(0.01.*wo)-1); %conversion to molality from wt%
            out = SeaFreeze([P_MPa' T_K' mo*ones(length(T_K),1)],'NH3');
            vel_kms = out.vel';
    end
end % fluidSoundSpeeds
function zero_alpha = alphaAdjust(P_MPa,T_K,wo,comp)
    [~,~,zero_alpha]= fluidEOS(P_MPa,T_K,wo,comp);
end % alphaAdjust
%% Core Size
function [C2MR2,R_fe] = CoreSize(rho_s,rho_fe,C_H2O,M_above_kg,R_s, M_kg,R_m)
    try
        R_fe = fzero(@(R_fe) getR_fe(rho_s,rho_fe,M_above_kg,R_s,R_fe, M_kg),[0 R_m]);
    catch
        R_fe = NaN;
    end
    C2MR2 = C_H2O+8/15*pi*((R_s.^5-R_fe.^5)*rho_s+rho_fe*R_fe.^5);%/Planet.M_kg/Planet.R_m^2;
end %CoreSize
function zero_me = getR_fe(rho_s,rho_fe,M_above_kg,R_s,R_fe, M_kg)
    zero_me = 4*pi/3*rho_fe*R_fe.^3 - (M_kg - M_above_kg-4*pi/3*rho_s*(R_s.^3-R_fe.^3));
end %getR_fe
function intout = interp1nan(x,y) % perplex outputs for VP, VS, KS, and GS seem to often have a few nans
    ninds = isnan(y);
    intout = interp1(x(~ninds),y(~ninds),x);
    if isnan(intout(end))
        intout(end)=intout(end-1);
    end
end %interp1nan
            
%% formatting for output
function d_str = getTableStr(Tb,Xin)
    d_str = {};
    for iT = 1:length(Tb)
        if Xin(iT)
            if mod(Xin(iT),1)
                d_str{iT} = [' &' num2str(Xin(iT),'%0.2f') ];
            else
                d_str{iT} = [' &' num2str(Xin(iT),'%0.0f') ];
            end
        else
            d_str{iT} = '& -';
        end
    end
end % getTableStr
function mantle = loadMantleEOS(str_meos)
[~,mantle]=read_perplex_table_nvars(str_meos);
end %loadMantleEOS
function mfluids = loadMantleFluids(str_mfluids)
[~,mfluids]=read_perplex_table_nvars(str_mfluids);
% [~,mfluids]=read_perplex_mantlefluids_table(str_mfluids);
%
% mfluids.co2_fn = griddedInterpolant(1e3*mfluids.p,mfluids.t,mfluids.co2);%
% mfluids.ch4_fn = griddedInterpolant(1e3*mfluids.p,mfluids.t,mfluids.ch4);
% mfluids.h2s_fn = griddedInterpolant(1e3*mfluids.p,mfluids.t,mfluids.h2s);
% mfluids.h2_fn = griddedInterpolant(1e3*mfluids.p,mfluids.t,mfluids.h2);
% mfluids.h2o_fn = griddedInterpolant(1e3*mfluids.p,mfluids.t,mfluids.h2o);
end %loadMantleFluids
function mphases = loadMantlePhases(str_mphases)
[~,mphases]=read_perplex_table_nvars(str_mphases);
end %loadMantlePhases

%% Plotting
function plotSolidInterior(prop,title_str,T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
if strcmp(opts.Punits,'GPa')
    xp = 1e-3;
else
    xp = 1;
end

    % correct for rotated p and t dimensions due to different versions of
    % Perple_X. Similar code is inserted a few lines down for the core. This change might break ungracefully if input silicate or core properties are
    % on a square matrix
    [m,n] = size(mantle.t);
    plotme = reshape(mantle.(prop),m,n);
if strcmp(prop,'Ks') || strcmp(prop,'Gs')
    pcolor(mantle.t,mantle.p*1e3*xp,plotme*1e-4);
else
    pcolor(mantle.t,mantle.p*1e3*xp,plotme);
end

hold on;
if ~isempty(core)
    [m,n] = size(core.t); % new code
    propcore = core.([prop '_fn'])(pcore,tcore); % new code
    if strcmp(prop,'Ks') || strcmp(prop,'Gs')
        pcolor(tcore,pcore*xp,propcore*1e-4);
    else
        pcolor(tcore,pcore*xp,propcore);
    end
else
    propcore = 0;
end
hp = plot(T_Planet_K',P_Planet_MPa'*xp,'w-','LineWidth',1);
maxcore = max(max(propcore));
maxmantle = max(max(plotme));
if maxcore>maxmantle
    maxprop = maxcore;
else
    maxprop = maxmantle;
end
if strcmp(prop,'Ks') || strcmp(prop,'Gs')
    maxprop = maxprop*1e-4;
elseif strcmp(prop,'alpha')
    maxprop = 5e-5;
end
shading interp; hc = colorbar; caxis([0 maxprop]);
box on; set(gca,'ydir','reverse');
if opts.Ttight
    xlims = get(gca,'XLim');
    set(gca,'XLim',[xlims(1) max(max(T_Planet_K))*1.01])
end
ylabel(['Pressure (' opts.Punits ')'],'Fontsize',14);
xlabel('Temperature (K)','Fontsize',14);
title(title_str,'Fontsize',18)

end %plotSolidInterior

function t_model=create_Taup(inputfile);
fid=fopen(inputfile);
line_num=1;
if contains(inputfile,'Titan')
    planet_Rad=2574.7;
end
rad_prev=0;
rho_prev=0;
Vp_prev=0;
Vs_prev=0;
phase_prev=0;
bmfile=strrep(inputfile,'.txt','.tvel')
fileID = fopen(bmfile,'w');
lay=1;
while ~feof(fid)
    tline = fgetl(fid);
    if line_num==1
         
        fprintf(fileID,'%s -P \n',inputfile)
        fprintf(fileID,'%s -S \n',inputfile)
   
    else
        contains_phase=contains(tline,'phase');
        if contains(tline,'NaN') % for last line in Planet Profile
            depth=planet_Rad; % read in radius and convert to m
            rho=rho; % uses density of previous line
            Vp=Vp; % grabs velocity and applies conversion
            Vs=Vs;
            phase=5;
        else
            depth=planet_Rad-str2num(tline(25:35)); % read in radius and convert to m
            rho=str2num(tline(37:47))/1000; % grab density
            Vp=str2num(tline(49:59)); % grabs velocity and applies conversion
            Vs=str2num(tline(61:71));
            Qmu=str2num(tline(73:83));
            L=4/3*(Vs/Vp).^2;
            Qp=(1/L)*Qmu;
            if L==0
                Qp=10000.0;
            end
        end
        
        
        if Vs==0 & lay==1; % change to Vs=0 indicating ocean
            lay=2;
            %Qmu=0;
            fprintf(fileID,'%.3f \t %.4f \t %.4f \t %.4f \n',...
                depth,Vp_prev,Vs_prev,rho_prev);
            % fprintf(fileID,'mantle \n')
        elseif Vs>0 && lay==2; % silicate core
            fprintf(fileID,'%.3f \t %.4f \t %.4f \t %.4f  \n',...
                depth,Vp_prev,Vs_prev,rho_prev);
            %fprintf(fileID,'outer-core \n')
            lay=3;
            cmb=line_num;
            %Qmu=1000;
        elseif (Qp>10000 | Qmu>10000) & lay==3
            fprintf(fileID,'%.3f \t %.4f \t %.4f \t %.4f \n',...
                depth,Vp_prev,Vs_prev,rho_prev);
            %fprintf(fileID,'inner-core \n')
            lay=4;
           
        end
        if (Vs~=0 | lay~=3) % issue with sometimes have an odd line where Vs is briefly non zero
        fprintf(fileID,'%.3f \t %.4f \t %.4f \t %.4f \n',...
            depth,Vp,Vs,rho);
        end
        depth_prev=depth;
        rho_prev=rho;
        Vp_prev=Vp;
        Vs_prev=Vs;
        Qp_prev=Qp;
        Qmu_prev=Qmu;
    end
    line_num=line_num+1;
end


fclose(fileID);
disp('TauP File Saved. Creating TauP Model. If this process takes more than a minute consider raising ice-ocean temperature or skipping this step')
try
    taup_model=taupcreate(bmfile);
catch
    disp('Adding small layer to fix issue with ray param')
    fid=fopen(bmfile);
    C=textscan(fid,'%s','delimiter','\n');
    line_1=char(C{1,1}(3));
    line_2=char(C{1,1}(4));
    [N M]=size(C{1,1});
    Cnew=C;
    Cnew{1,1}(5:N+1)=C{1,1}(4:end);
    new_line=[line_2(1:5) line_1(6:end)];
    fileID = fopen(bmfile,'w');
    for k=1:numel(Cnew{1,1})
        if k==4
            fprintf(fileID,'%s \n',new_line);
        else
            fprintf(fileID,'%s \n',Cnew{1,1}{k});
        end
    end
    fclose(fileID);
    try
    taup_model=taupcreate(bmfile);
    catch
        disp('Error creating Taup Model. Try adding more layers in the near-surface or manually adjusting velocity profile')
    end
    taupcurve('mod',taup_model)
end
t_model=taup_model;
end