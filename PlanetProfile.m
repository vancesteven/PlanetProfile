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
        
if strcmp(Planet.Ocean.comp,'Seawater')
    wtPpt = Planet.Ocean.w_ocean_pct; % We are already using ppt for Seawater even though the variable is called WtPct
else
    wtPpt = 10*Planet.Ocean.w_ocean_pct;
end
       
    %% steve to clean this up by calling Params.cfg
Params.NOPLOTS = cfg.NO_PLOTS;
Params.CALC_NEW = cfg.CALC_NEW;
Params.CALC_NEW_REFPROFILES = cfg.CALC_NEW_REF;
Params.CALC_NEW_SOUNDSPEEDS = cfg.CALC_NEW_SOUND;
Params.INCLUDE_ELECTRICAL_CONDUCTIVITY = cfg.CONDUCT;
Params.HOLD = cfg.HOLD;
Params.savefigformat = cfg.fig_fmt;

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
[n_clath, n_iceI, n_ocean] = deal(zeros(1));
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
[z_m,r_m,g_ms2,M_above_kg,M_below_kg] = deal(zeros(1,nsteps));
M_above_kg(1) = 0;
M_below_kg(1) = Planet.M_kg;
r_m(1) = Planet.R_m;
g_ms2(1) = Gg*Planet.M_kg/Planet.R_m^2;
D_conductivityIh = 632; % W m-1; Andersson et al. 2005 (For comparison, Mckinnon 2006 uses a value of 621 from Slack 1980)

[Planet.Profile_fname, Planet.Profile_ID] = deal(strings(1));

% Preallocate for ConvectionDeschampsSotin
[Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I,Qb,Planet.zb_outerIce_m] = deal(zeros(1));


%%
%if ~Planet.NoH2O
    if Params.CALC_NEW
        % adds number of clathrates based on inputs or sets to zero, also sets
        % maximum depth if specified

%%      

        [T_K,P_MPa,rho_kgm3, phi_hydro_frac] = deal(zeros(1,nsteps));

        phase = zeros(1,nsteps);
        % assign phase based on clathrates or Ih
        phase((1+Params.nsteps_clath):(Params.nsteps_clath+Params.nsteps_iceI))=1;
        phase(1:Params.nsteps_clath)=30; % index to indicate clathrates, which are not an ice phase
        Tb_K = Planet.Tb_K;
    

        %% IceLayers --- function iteratively setting up the thermal profile, the density and temperature of the layer with each pressure step
        %ice Ih, ice III, ice V
        %--------------------------------------------------------------------------

        T_K(1) = Planet.Tsurf_K;
        P_MPa(1) = Planet.Psurf_MPa;
        % if ice shell had to be thinned in previous run, this will
        % reset indexing correctly. Set starting values:
        n_iceI = Params.nsteps_iceI; % these are redundant. just change text below to use Params.n...
        n_clath = Params.nsteps_clath;
        n_ocean = Params.nsteps_ocean; % 

        if n_clath>0
            clath_out{1} = Helgerud_sI(P_MPa(1),T_K(1));
            rho_kgm3(1)=clath_out{1}.rho;
            [Cp(1) alpha_K(1)]= getCpIce(P_MPa(1),T_K(1),phase(1)) ;
        else
            rho_kgm3(1) = getRhoIce(P_MPa(1),T_K(1),1);
            [Cp(1) alpha_K(1)]= getCpIce(P_MPa(1),T_K(1),phase(1)) ;
        end
        try
            Pb_MPa = getPfreeze(Planet.Tb_K,wo,Planet.Ocean.comp);
            %[Cp(1) alpha_K(1)]= getCpIce(P_MPa(1),T_K(1),phase(1)) ;


            deltaP = Pb_MPa/(n_clath+n_iceI-1);
            % assumes Pressure gradient doesnt change betweeen clathrates
            % and need to check assumption
            % if n_clath == 0 this will skip and move on to next section
            %% ClathrateLayer function within IceLayers
            for il=2:n_clath % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                P_MPa(il) = P_MPa(il-1) + deltaP;
                T_K(il) = (Planet.Tb_K.^(P_MPa(il)./Pb_MPa)).*(Planet.Tsurf_K.^(1-P_MPa(il)./Pb_MPa));
                clath_out{il} = Helgerud_sI(P_MPa(il),T_K(il));

                rho_kgm3(il)=clath_out{il}.rho;
                [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),phase(il)) ;

                z_m(il) = z_m(il-1)+ (P_MPa(il)-P_MPa(il-1))*1e6/g_ms2(il-1)/rho_kgm3(il-1);
                r_m(il) = Planet.R_m-z_m(il);

                % determine local gravity and check to make sure maximum
                % clathrate depth isn't reached
                M_below_kg(il) = M_above_kg(il-1) + 4/3*pi*(r_m(il-1)^3-r_m(il)^3)*rho_kgm3(il);
                M_below_kg(il) = Planet.M_kg-M_above_kg(il);
                g_ms2(il) = Gg*M_below_kg(il)/r_m(il)^2;
                if z_m(il)>max_clath_depth % check to see if maximum depth of clathrates reached
                    n_iceI=n_iceI+(n_clath-il)+1;
                    n_clath=il-1;
                    phase((1+n_clath):(n_clath+n_iceI))=1;
                    break
                end
                %rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1);
            end % end ClathrateLayer
            if n_clath>0
               ii=n_clath+1;
            else
                ii=2;
            end

            %% IceILayer
            for il=ii:(n_iceI+n_clath) % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                P_MPa(il) = P_MPa(il-1) + deltaP;
                T_K(il) = (Planet.Tb_K.^(P_MPa(il)./Pb_MPa)).*(Planet.Tsurf_K.^(1-P_MPa(il)./Pb_MPa));
                rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1); % I THINK THIS CALL CAN ACCEPT A VECTOR
                [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),phase(il)) ; % I THINK THIS CALL CAN ACCEPT A VECTOR

                %rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1);
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
                phase((n_iceI+n_clath)-nIceIIILithosphere:(n_iceI+n_clath))=3;

                PbI_MPa = 210; % the Ih-III transition is essentially fixed, within a few MPa
                deltaP = PbI_MPa/((n_iceI+n_clath)-5-1);  % save five entries at the bottom for ice III

                for il=2:n_clath % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                P_MPa(il) = P_MPa(il-1) + deltaP;
                T_K(il) = (Planet.Tb_K.^(P_MPa(il)./Pb_MPa)).*(Planet.Tsurf_K.^(1-P_MPa(il)./Pb_MPa));
                clath_out{il} = Helgerud_sI(P_MPa(il),T_K(il));

                rho_kgm3(il)=clath_out{il}.rho;
                 z_m(il) = z_m(il-1)+ (P_MPa(il)-P_MPa(il-1))*1e6/g_ms2(il-1)/rho_kgm3(il-1);
                 r_m(il) = Planet.R_m-z_m(il);
                 [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),phase(il))

                 % determine local gravity
                 M_above_kg(il) = M_above_kg(il-1) + 4/3*pi*(r_m(il-1)^3-r_m(il)^3)*rho_kgm3(il);
                 M_below_kg(il) = Planet.M_kg-M_above_kg(il);
                 g_ms2(il) = Gg*M_below_kg(il)/r_m(il)^2;
                 if z_m(il)>max_clath_depth % check to see if maximum depth of clathrates reached
                    n_iceI=n_iceI+(n_clath-il)+1;
                    n_clath=il-1;
                    phase(:,(1+n_clath):(n_clath+n_iceI))=1;
                    break
                 end 

                %rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1);
                end
                if n_clath>0
                    ii=n_clath+1;
                else
                    ii=2;
                end

                for il=ii:(n_iceI+n_clath)-nIceIIILithosphere % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                   P_MPa(il) = P_MPa(il-1) + deltaP;
                 T_K(il) = (Planet.Tb_K.^(P_MPa(il)./Pb_MPa)).*(Planet.Tsurf_K.^(1-P_MPa(il)./Pb_MPa));
                    rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1);
                    [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),phase(il)) ;

                %rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1);
                end


                TbIII = Planet.Tb_K+1;
                Pb_MPa = getPfreezeIII(TbIII,wo,Planet.Ocean.comp);% fix the thickness of the ice III layer based on T>Tb by 2K
                deltaP = (Pb_MPa-PbI_MPa)/nIceIIILithosphere;  %  five entries at the bottom for ice III
                for il = (n_clath+n_iceI)-nIceIIILithosphere+1:(n_iceI+n_clath)
                    P_MPa(il) = P_MPa(il-1) + deltaP;
                    T_K(il) = (TbIII.^(P_MPa(il)./Pb_MPa)).*(TbIII.^(1-P_MPa(il)./Pb_MPa));
                    rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),3);
                    [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),3) ;
                end

            elseif isfield(Params,'BOTTOM_ICEV') && Params.BOTTOM_ICEV
                disp('Adding ice V and III to the bottom of ice Ih. Make sure the ocean salinity is high enough that doing this makes sense')
                nIceIIILithosphere=5;
                nIceVLithosphere=5;
                phase(n_iceI-nIceVLithosphere-nIceIIILithosphere:n_iceI-nIceVLithosphere)=3;
                phase(n_iceI-nIceVLithosphere:n_iceI)=5;

                PbI_MPa = 210; % the Ih-V transition is essentially fixed, within a few MPa
                deltaP = PbI_MPa/((n_iceI+n_clath)-5-1);  % save five entries at the bottom for ice V

                %clathrates
                for il=2:n_clath % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                      P_MPa(il) = P_MPa(il-1) + deltaP;
                      T_K(il) = (Planet.Tb_K.^(P_MPa(il)./Pb_MPa)).*(Planet.Tsurf_K.^(1-P_MPa(il)./Pb_MPa));
                      clath_out{il} = Helgerud_sI(P_MPa(il),T_K(il));

                      rho_kgm3(il)=clath_out{il}.rho;
                      [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),phase(il))
                      z_m(il) = z_m(il-1)+ (P_MPa(il)-P_MPa(il-1))*1e6/g_ms2(il-1)/rho_kgm3(il-1);
                      r_m(il) = Planet.R_m-z_m(il);

                      % determine local gravity
                      M_above_kg(il) = M_above_kg(il-1) + 4/3*pi*(r_m(il-1)^3-r_m(il)^3)*rho_kgm3(il);
                      M_below_kg(il) = Planet.M_kg-M_above_kg(il);
                      g_ms2(il) = Gg*M_below_kg(il)/r_m(il)^2;
                      if z_m(il)>max_clath_depth % check to see if maximum depth of clathrates reached
                          n_iceI=n_iceI+(n_clath-il)+1;
                          n_clath=il-1;
                          phase(:,(1+n_clath):(n_clath+n_iceI))=1;
                          break
                    end

                    %rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1);
                end
                if n_clath>0
                    ii=n_clath+1;
                else
                    ii=2;
                end

                %ice Ih
                % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                for il=ii:(n_iceI+n_clath)-nIceVLithosphere % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                    P_MPa(il) = P_MPa(il-1) + deltaP;
                    T_K(il) = (Planet.Tb_K.^(P_MPa(il)./Pb_MPa)).*(Planet.Tsurf_K.^(1-P_MPa(il)./Pb_MPa));
                    rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),phase(ill));
                    [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),phase(ill)) ;
                    %rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),1);
                end

                %ice III
                TbIII = Planet.Tb_K+1;
                Pb_MPa = getPfreezeIII(TbIII,wo,Planet.Ocean.comp);% fix the thickness of the ice III layer based on T>Tb by 2K
                deltaP = (Pb_MPa-PbI_MPa)/(nIceIIILithosphere);  %  five entries at the bottom for ice III
                for il = (n_iceI+n_clath)-nIceIIILithosphere+1:(n_iceI+n_clath)
                    P_MPa(il) = P_MPa(il-1) + deltaP;
                    T_K(il) = (TbIII.^(P_MPa(il)./Pb_MPa)).*(TbIII.^(1-P_MPa(il)./Pb_MPa));
                    rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),3);
                    [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),3) ;
                end

                %ice V
                TbV = Planet.Tb_K+1;
                Pb_MPa = getPfreezeV(TbV,wo,Planet.Ocean.comp);% fix the thickness of the ice V layer based on T>Tb by 2K
                deltaP = (Pb_MPa-PbI_MPa)/(nIceVLithosphere);  %  five entries at the bottom for ice V
                for il = (n_iceI+n_clath)-nIceVLithosphere+1:(n_iceI+n_clath)
                    P_MPa(il) = P_MPa(il-1) + deltaP;
                    T_K(il) = (TbV.^(P_MPa(il)./Pb_MPa)).*(TbV.^(1-P_MPa(il)./Pb_MPa));
                    rho_kgm3(il) = getRhoIce(P_MPa(il),T_K(il),3);
                    [Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),3) ;
                end

                %adjusts for surface porosity if set
                if POROUS_ICE % Steve to ask Angela if this is actually being used.
                    % correction for porosity
                    por_in.p=P_MPa(1:il)*1e-3;
                    por_in.t = T_K(1:il);
                    por_in.den = rho_kgm3(1:il);
                    %por_in.vp = velsIce.Vclathl_kms;
                    %por_in.vs = velsIce.Vclatht_kms;
                    if isfield(Planet,'phi_surface')
                        por_out = get_porosity_ice(por_in,Planet.phi_surface);
                    else
                        por_out =get_porosity_ice(por_in);
                    end
                    %permeability = por_out.per;
                    %ice_ind=find(phase==30);
                    rho_kgm3(1:il) = por_out.den;
                    phi_hydro_frac(1:il)=por_out.por;
                    %velsIce.Vclathl_kms = por_out.vp;
                    %velsIce.Vclatht_kms = por_out.vs;

                    disp(['Average ice porosity: ' num2str(mean(phi_hydro_frac(1:il)))])
                    disp(['Porosity: ' num2str(phi_hydro_frac(1:il))])
                end
            end % end IceIIIUnderplateLayer
        end

            %% OceanLayer
        %OCEAN + ICE III/V/VI SHELL
        %--------------------------------------------------------------------------
        deltaP = (Params.Pseafloor_MPa-Pb_MPa)/n_ocean; %
        if deltaP<=0 % prevent an embarassing error that can occur when adapting a file for a small object to that of a larger one.
            error('negative increment of pressure while computing ocean profile. Be sure Params.Pseafloor_MPa is greater than the likely pressure at the silicate interface. Currently it''s less than the pressure at the base of the ice I layer.')
        end
        for il =1:n_ocean
            ill = il+n_iceI+n_clath;
            disp(['il: ' num2str(il) '; P_MPa: ' num2str(round(P_MPa(ill-1))) '; T_K: ' num2str(round(T_K(ill-1)))]);
            P_MPa(ill) = Pb_MPa + il*deltaP;

            if il==1 % establish the phase vector
                phase(ill) = getIcePhase(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
            elseif phase(ill-1)~=6
                phase(ill) = getIcePhase(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
                if phase(ill-1)==3 && phase(ill)==0 % fix to instabilities in the phase lookup for ice III
                    disp('Fixing apparent instability in the EOS?ice III is forming above the ocean where it wasn''t specified')
                    phase(ill)=3;
                end
                if phase(ill-1)==0 && phase(ill)==1 % fix to instabilities in the phase lookup for ice I
                    disp('Fixing apparent instability in the EOS?ice Ih is forming under the ocean')
                    phase(ill)=0;
                end
                if phase(ill-1)==5 && phase(ill)==0 % fix to instabilities in the phase lookup for ice V for Titan NH3 LBF svance Jul 6 2021
                    disp('Fixing apparent instability in the EOS? fluid is forming under ice V')
                    phase(ill)=5;
                end
            else
                phase(ill) = 6;
            end

            if phase(ill) == 0 %if ocean
                [rho_ocean,Cp(ill),alpha_o]= fluidEOS(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);
                % Forbidding MgSO4 in the check below avoids a problem in
                % application of the EOS for low-salinity MgSO4 oceans.
                % Negative thermal expansivity regions in the MgSO4 EOS may be
                % artifacts of the current EOS calculation. See Vance et al. 2014
                if alpha_o<=0 && ~strcmp(Planet.Ocean.comp,'MgSO4') && ~Planet.ALLOW_NEGALPHA
                    disp('Ocean alpha at ice interface is less than zero. adjusting temperature upward.')
                    disp('This means there''s a conductive layer at the interface with thickness inversely proportional to the heat flow.')
                    disp('The thickness is likely less than a few 100 m. See Melosh et al. 2004.')
                    Tnew = fzero(@(T) alphaAdjust(P_MPa(ill),T,wo,Planet.Ocean.comp),273);
                    disp(['The new temperature is ' num2str(Tnew) ', or delta T = ' num2str(Tnew-T_K(ill-1)) '.'])
                    T_K(ill)=Tnew;
                else
                    T_K(ill) = T_K(ill-1)+ alpha_o*T_K(ill-1)./(Cp(ill))./(rho_ocean)*deltaP*1e6; %adiabatic gradient in ocean; this introduces an error, in that we are using the temperature from the previous step
                end
                alpha_K(ill) = alpha_o;
                rho_kgm3(ill) = fluidEOS(P_MPa(ill),T_K(ill),wo,Planet.Ocean.comp); % to be deleted

%                         [rho_ocean,Cp(ill),alpha_o]= fluidEOS(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);% THIS IS REDUNDANT TO THE CALCULATION OF RHO_OCEAN ABOVE
            else %% HPIceLayers
                % allow for stable dense fluid under high pressure ice --> this
                % was added for the callisto study. Commented out currently
                rhoice = getRhoIce(P_MPa(ill),T_K(ill-1),phase(ill));
                % rhoocean = fluidEOS(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);
                % if rhoocean>=rhoice
                %     phase(ill)=0;
                % end
                if ~isfield(Planet.Ocean,'fnTfreeze_K')
                    if il==6
                        x = 1;
                    end
                    T_K(ill) = getTfreeze(P_MPa(ill),wo,Planet.Ocean.comp,T_K(ill-1)); %find the temperature on the liquidus line corresponding to P(k,i+nIceI); should really use a conductive profile here, but that would seem to naturally bring us back to the liquidus. Steve wants to confirm.
                else
                    T_K(ill) = Planet.Ocean.fnTfreeze_K(P_MPa(ill),wo);
                end
                if T_K(ill)<T_K(ill-1)
                    T_K(ill)=T_K(ill-1); % this may no longer be needed because negalpha is accounted for above -- sincerely, steve 6/15/21
                end
                rho_kgm3(ill) = getRhoIce(P_MPa(ill),T_K(ill),phase(ill));

%                     [rho_ocean,Cp(ill),alpha_o]=
%                     fluidEOS(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);
%                     % this is odd and doesn't seem to be needed--steve
%                     6/15/21
                %[Cp(il) alpha_K(il)]= getCpIce(P_MPa(il),T_K(il),phase(ill)) ;
            end
        end

        % MJS 2021-10-29: I don't think these lines are important or do
        % anything substantial. Both quantities are recalculated
        % aplenty in the next loops. This was being done just before
        % saving/reloading in the master branch before now.
        rho_kgm3(1) = rho_kgm3(2); % continuity
        Cp(1)=Cp(2);


%%%%%%%%%%%%%%%%%%%%%
% convert to depth —- PlanetDepths
%%%%%%%%%%%%%%%%%%%
%% calculate gravity in each layer instead of assuming surface gravity applies.
% allocate variables
%% HydrosphereDepths
    deltaP = Pb_MPa/(n_iceI+n_clath);

    % calculates depth for clathrates and ice separatly so the depths
    % can be recorded accurately
    for il = 2:(n_clath)
        % calculate depth
        z_m(il) = z_m(il-1)+ (P_MPa(il)-P_MPa(il-1))*1e6/g_ms2(il-1)/rho_kgm3(il-1);
        % convert to radius
        r_m(il) = Planet.R_m-z_m(il); 

        % determine local gravity
        M_above_kg(il) = M_above_kg(il-1) + 4/3*pi*(r_m(il-1)^3-r_m(il)^3)*rho_kgm3(il);
        M_below_kg(il) = Planet.M_kg-M_above_kg(il);
        g_ms2(il) = Gg*M_below_kg(il)/r_m(il)^2;
    end
    %checks if there were clathrates or not
    if n_clath>0
        ii=n_clath+1;
        Zclath=z_m(n_clath);
        Planet.Zclath=Zclath;
    else
        ii=2;
        Zclath=0;
        Planet.Zclath=0;
    end
    for il = ii:(n_clath+n_iceI)
        % calculate depth
        z_m(il) = z_m(il-1)+ (P_MPa(il)-P_MPa(il-1))*1e6/g_ms2(il-1)/rho_kgm3(il-1);
        % convert to radius
        r_m(il) = Planet.R_m-z_m(il); 

        % determine local gravity
        M_above_kg(il) = M_above_kg(il-1) + 4/3*pi*(r_m(il-1)^3-r_m(il)^3)*rho_kgm3(il);
        M_below_kg(il) = Planet.M_kg-M_above_kg(il);
        g_ms2(il) = Gg*M_below_kg(il)/r_m(il)^2;
    end
    if isempty(il)
        zb_outerIce_m=z_m(ii-1);
        Planet.zb_outerIce_m=zb_outerIce_m;
    else
        zb_outerIce_m=z_m(il);
        Planet.zb_outerIce_m=zb_outerIce_m;
    end
    deltaP = (Params.Pseafloor_MPa-Pb_MPa)/n_ocean; %
    for il = 1:n_ocean
        ill = il+(n_iceI+n_clath);
        %calculate depth
        dz = deltaP*1e6/g_ms2(il-1+(n_iceI+n_clath))/rho_kgm3(ill);
        z_m(ill) = z_m(il-1+(n_iceI+n_clath))+ dz; % using the previous gravity step, since we haven't calculated the present step.  this introduces an error
        % convert to radius
        r_m(ill) = Planet.R_m-z_m(ill); 

        % determine local gravity
        M_above_kg(ill) = M_above_kg(il-1+(n_iceI+n_clath))+4/3*pi*((r_m(ill)+dz)^3-r_m(ill)^3)*rho_kgm3(ill);
        M_below_kg(ill) = Planet.M_kg-M_above_kg(ill);
        g_ms2(ill) = Gg*M_below_kg(ill)/r_m(ill)^2;
    end
    disp(['z_iceI: ' num2str(zb_outerIce_m/1e3) ' km'])
    z_ocean_m= z_m(n_ocean); % depth to the ocean



    %% compute conductive heat through the ice I layer
    Qb = D_conductivityIh*log(Planet.Tb_K/Planet.Tsurf_K)/Planet.zb_outerIce_m;

    %% compute solid state convection ice
    % We use these values much later, but we need them to construct
    % the filename shortly, so we calculate them now. This allows
    % us to skip unnecessary calculations in the case of
    % cfg.SKIP_PROFILES=1
      if max_clath_depth<Planet.zb_outerIce_m % only a clathrate lid
        % asummes Q across ice-ocean is same Q across clathrates/ice. 
       [Q_Wm2, T_clath_ice, deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]= ...
   clathrate_lid_thermo(Planet.Tsurf_K,Planet.Tb_K,P_MPa(1:n_clath+n_iceI), n_clath,n_iceI,Planet.zb_outerIce_m, max_clath_depth,g_ms2(1));      
      else
       [Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]=...
        ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K,PbI_MPa/2,Planet.zb_outerIce_m,g_ms2(1),phase(1)+1);
     end

    % Make this calculation now in order to get Planet.Qmantle_Wm2 for making
    % filenames shortly       
    if CONVECTION_FLAG_I && eTBL_m>zb_outerIce_m
        warning('Convection predicted by not possible becuase the conductive layer thickness exceeds the thickness of the ice.')
        disp('Perhaps T_surf is outside the valid range for the scaling from Deschamps and Sotin 2001.')
        disp('Setting CONVECTION_FLAG_I to zero')
        CONVECTION_FLAG_I = 0;
    end
    % Convection has to be calculated prior to assigning depths in case
    % ice shell needs to be thinned to account for clathrates

    if CONVECTION_FLAG_I
        %conductive upper layer
        nConvectIce=(n_iceI+n_clath)-nIceIIILithosphere-1; % indices where ice/claths exist in upper layer
        nconold=nConvectIce;
        if max_clath_depth<Planet.zb_outerIce_m
            %eTBL_m=eTBL_m+max_clath_depth; fixed in new
            %version
            T_K(1:n_clath)=linspace(Planet.Tsurf_K,T_clath_ice,n_clath);
            inds_eTBL = find(z_m(n_clath:nConvectIce)<=eTBL_m);

            if ~isempty(inds_eTBL)
                inds_eTBL=inds_eTBL +n_clath;
                if length(inds_eTBL)>1
            Pterm =(P_MPa(inds_eTBL)-P_MPa(inds_eTBL(1)))./(P_MPa(inds_eTBL(end))-P_MPa(inds_eTBL(1)));

            T_K(inds_eTBL) = (Tc.^(Pterm)).*(T_clath_ice.^(1-Pterm));
                else
                     T_K(inds_eTBL)=T_clath_ice+Q_Wm2./kIce.*(z_m(inds_eTBL)-z_m(inds_eTBL-1));
                end
            else
                inds_eTBL=n_clath;
            end
        else
        inds_eTBL = find(z_m(1:nConvectIce)<=eTBL_m);

        Pterm = P_MPa(inds_eTBL)./P_MPa(inds_eTBL(end));
        if phase(inds_eTBL(1))==1
        T_K(inds_eTBL) = (Tc.^(Pterm)).*(Planet.Tsurf_K.^(1-Pterm));
        elseif phase(inds_eTBL(1))==30;
             T_K(1:inds_eTBL(end))=linspace(Planet.Tsurf_K,Tc,length(inds_eTBL));
        end
        end

         rho_kgm3(1:inds_eTBL(end))=getRhoIce(P_MPa(1:inds_eTBL(end)),T_K(1:inds_eTBL(end)),phase(1:inds_eTBL(end)));

        P_bound=P_MPa(inds_eTBL(end));
        T_K_ccbound=T_K(inds_eTBL(end));


        %convective region

        for iconv = inds_eTBL(end)+1:nConvectIce % added parentheses around 1:nConvectIce 20200103 SDV

            rho_kgm3(iconv)=getRhoIce(P_MPa(iconv),T_K(iconv-1),phase(iconv));

            if POROUS_ICE % adjust if porosity needs to be considered
                por_in.p=P_MPa(iconv)*1e-3;
                por_in.t = T_K(iconv-1);
                por_in.den = rho;

                if isfield(Planet,'phi_surface')
                    por_out = get_porosity_ice(por_in,Planet.phi_surface);
                else
                    por_out =get_porosity_ice(por_in);
                end
                rho = por_out.den;
                phi = pot_out.por;
            else
                phi = 0;
            end
                %[Cp(iconv), alpha_K(iconv)] = getCpIce(P_MPa(iconv),T_K(iconv-1),1);
                if phase(iconv)==1 % if water ice, use SeaFreeze otherwise use Helgeurd.
                    cpout=SeaFreeze([P_MPa(iconv),T_K(iconv-1)],'Ih');
                    alpha_K(iconv)=cpout.alpha;
                    Cp(iconv)=cpout.Cp;
                else

                    cpout=Helgerud_sI(P_MPa(iconv),T_K(iconv-1));
                    [Cp(iconv) alpha_K(iconv)]= getCpIce(P_MPa(iconv),T_K(iconv-1),phase(iconv)) ;
                    alpha_K(iconv)=cpout.alpha; % use better alpha
                end

            %aK = 1.56e-4; % thermal expansive? Switch and use SeaFreeze ( find clath values()

            T_K(iconv) = T_K(iconv-1)+alpha_K(iconv).*T_K(iconv)./Cp(iconv)./rho_kgm3(iconv)*deltaP*1e6;
            % double check temperatures make sense
            if strcmp(Planet.Ocean.comp,'NH3') % kluge svance july 7 2021. not sure why two phase tests are happening. generally, the seafreeze test should be performed when possible
                phase_test = 1;
            else
                phase_test = getIcePhase(P_MPa(iconv),T_K(iconv-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
            end
            phase_test2=SF_WhichPhase([P_MPa(iconv),T_K(iconv)]);

            if phase_test==0 || phase_test2==0 % & n_clath>0
                Planet.Tb_K=T_K(iconv);
                %iconv=iconv-2;
                Zm_break=z_m(iconv); % saves the depth at which is broked
                zb_outerIce_m=Zm_break;   % new ice depth is saved

                Zdiff=z_m(nConvectIce)-Zm_break; % change in depth
                 % previous convection layer

                nConvectIce=iconv;% changes the indices for convection


                %make adjustments to indices n_clath+n_ice+n_ocean
                %should always be the same.
                n_iceI = iconv - n_clath; % new ice layer is the last index-n_clath
                if n_iceI<0
                    n_iceI = 0; % indicates entire ice shell would be clathrates
                    iconv = iconv-1;
                    Zclath = z_m(iconv);
                end

                n_ocean = nsteps - n_iceI - n_clath;% adds more indices to ocean


                PbI_MPa=P_MPa(iconv);
                Pb_MPa=P_MPa(iconv);
                %                        alculate parameters
%                     [Q_Wm2_new,deltaTBL_m_new,eTBL_m_new,Tc_new,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]=...
%                         ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K,PbI_MPa/2,zb_outerIce_m,g_ms2(1),phase(iconv)+1);
                inds_deltaTBL = find(z_m(1:nConvectIce)>=z_m(nConvectIce)-deltaTBL_m);
                T_K(inds_deltaTBL) = (Planet.Tb_K.^(P_MPa(inds_deltaTBL)./PbI_MPa)).*(T_K(inds_deltaTBL(1)-1).^(1-P_MPa(inds_deltaTBL)./PbI_MPa));

                % %
                %
                deltaP = (Params.Pseafloor_MPa-Pb_MPa)/n_ocean; % in case PbMPa changed
                % recalculate ocean
                for il=1:n_ocean
                    ill=il+inds_deltaTBL(end);
                    % adjust for new ocean layers if necessary
                    P_MPa(ill) = Pb_MPa + il*deltaP;
                    [rho_ocean,Cp(ill),alpha_o]= fluidEOS(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);
                    if alpha_o<=0

                        disp('Ocean alpha at ice interface is less than zero. adjusting temperature upward.')
                        disp('This means there''s a conductive layer at the interface with thickness inversely proportional to the heat flow.')
                        disp('The thickness is likely less than a few 100 m. See Melosh et al. 2004.')
                        Tnew = fzero(@(T) alphaAdjust(P_MPa(ill),T,wo,Planet.Ocean.comp),273);
                        disp(['The new temperature is ' num2str(Tnew) ', or delta T = ' num2str(Tnew-T_K(ill-1)) '.'])
                        T_K(ill)=Tnew;
                    else
                        T_K(ill) = T_K(ill-1)+ alpha_o*T_K(ill-1)./(Cp(ill))./(rho_ocean)*deltaP*1e6; %adiabatic gradient in ocean; this introduces an error, in that we are using the temperature from the previous step
                    end


                    phase(ill) = getIcePhase(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI

                    if phase(ill)==0
                        rho_kgm3(ill) = rho_ocean;
                         [rho_ocean,Cp(ill),alpha_o]= fluidEOS(P_MPa(ill),T_K(ill-1),wo,Planet.Ocean.comp);
                         alpha_K(ill) = alpha_o;

                    else
                        rho_kgm3(ill) = getRhoIce(P_MPa(ill),T_K(ill),phase(ill));

                        if ~isfield(Planet.Ocean,'fnTfreeze_K')

                            T_K(ill) = getTfreeze(P_MPa(ill),wo,Planet.Ocean.comp,T_K(ill-1)); %find the temperature on the liquidus line corresponding to P(k,i+nIceI); should really use a conductive profile here, but that would seem to naturally bring us back to the liquidus. Steve wants to confirm.
                        else
                            T_K(ill) = Planet.Ocean.fnTfreeze_K(P_MPa(ill),wo);
                        end
                        if T_K(ill)<T_K(ill-1)
                            T_K(ill) = T_K(ill-1);
                        end
                       [Cp(ill), alpha_K(ill)]= getCpIce(P_MPa(iconv),T_K(iconv-1),phase(ill)) ;


                    end

                end
                z_ocean_m=z_m(ill);

                break
            end
            phi_hydro_frac(iconv) = phi;
        end

        if nconold==nConvectIce
            % bottom layer of conduction, recalculates for phase at bottom
             [Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]=...
                    ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K,PbI_MPa/2,zb_outerIce_m,g_ms2(1),phase(iconv)+1);

            inds_deltaTBL = find(z_m(1:nConvectIce)>=z_m(nConvectIce)-deltaTBL_m);
            %
            %
            T_K(inds_deltaTBL) = (Planet.Tb_K.^(P_MPa(inds_deltaTBL)./PbI_MPa)).*(T_K(inds_deltaTBL(1)-1).^(1-P_MPa(inds_deltaTBL)./PbI_MPa));
            %           end
            rho_kgm3(inds_deltaTBL) = getRhoIce(P_MPa(inds_deltaTBL),T_K(inds_deltaTBL),phase(inds_deltaTBL));
             %[Cp(inds_deltaTBL) alpha_K(inds_deltaTBL)]= getCpIce(P_MPa(iconv),T_K(iconv-1),phase(inds_deltaTBL)) ;


            z_ocean_m= z_m(inds_deltaTBL(end)+1);

            if POROUS_ICE
                por_in.p=P_MPa(inds_deltaTBL)*1e-3;
                por_in.t = T_K(inds_deltaTBL);
                por_in.den = rho_kgm3(inds_deltaTBL);
                if isfield(Planet,'phi_surface')
                    por_out = get_porosity_ice(por_in,Planet.phi_surface);
                else
                    por_out =get_porosity_ice(por_in);
                end

                rho_kgm3(inds_deltaTBL) = por_out.den;
                phi_hydro_frac(inds_deltaTBL) = por_out.por;
            end
        end
    else
        z_ocean_m = zb_outerIce_m;

%             if nconold==nConvectIce
%                 % bottom layer of conduction, recalculates for phase at bottom
%                  [Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]=...
%                         ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K,PbI_MPa/2,zb_outerIce_m,g_ms2(1),phase(iconv)+1);
% 
%                 inds_deltaTBL = find(z_m(1:nConvectIce)>=z_m(nConvectIce)-deltaTBL_m);
%                 %
%                 %
%                 T_K(inds_deltaTBL) = (Planet.Tb_K.^(P_MPa(inds_deltaTBL)./PbI_MPa)).*(T_K(inds_deltaTBL(1)-1).^(1-P_MPa(inds_deltaTBL)./PbI_MPa));
%                 %           end
%                 rho_kgm3(inds_deltaTBL) = getRhoIce(P_MPa(inds_deltaTBL),T_K(inds_deltaTBL),phase(inds_deltaTBL));
%                  %[Cp(inds_deltaTBL) alpha_K(inds_deltaTBL)]= getCpIce(P_MPa(iconv),T_K(iconv-1),phase(inds_deltaTBL)) ;
%                         
% 
%                 z_ocean_m= z_m(inds_deltaTBL(end)+1);
% 
%                 if POROUS_ICE
%                     por_in.p=P_MPa(inds_deltaTBL)*1e-3;
%                     por_in.t = T_K(inds_deltaTBL);
%                     por_in.den = rho_kgm3(inds_deltaTBL);
%                     if isfield(Planet,'phi_surface')
%                         por_out = get_porosity_ice(por_in,Planet.phi_surface);
%                     else
%                         por_out =get_porosity_ice(por_in);
%                     end
% 
%                     rho_kgm3(inds_deltaTBL) = por_out.den;
%                 end
%             end

    end % if CONVECTION_FLAG


    if Planet.EQUIL_Q
        if (~isfield(Params,'NO_ICEI_CONVECTION') || Params.NO_ICEI_CONVECTION == false) && ...
           (CONVECTION_FLAG_I || (isfield(Params,'FORCE_ICEI_CONVECTION') && Params.FORCE_ICEI_CONVECTION == true))
            Qmantle_Wm2 = Q_Wm2;
        else
            Qmantle_Wm2 = Qb;
        end
    else
        Qmantle_Wm2 = 2.2e11/4/pi/Planet.R_m^2; % this is more reasonable for radiogenic only
    end
    Planet.Qmantle_Wm2 = Qmantle_Wm2;

    deltaP = (Params.Pseafloor_MPa-Pb_MPa)/n_ocean; %

    Planet.zClath_m = Zclath;

    Planet.Profile_ID = ['Ts' num2str(Planet.Tsurf_K,'%0.0f') 'Zb' strLow num2str(Planet.zb_outerIce_m,'%0.0f') ...
        'mQm' num2str(1000*Planet.Qmantle_Wm2,'%0.0f') 'mWm2_CMR2p' ...
        num2str(10000*Planet.Cmeasured,'%0.0f') '_' thiseos];
    Planet.Profile_fname = [savefile '_' char(Planet.Profile_ID)];
% end % ~Planet.NoH2O




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
npre = nz-nR; % confine the search for structure to the lower part % MJS 2021-10-31: This quantity does not make sense in how it is defined and used. It seems superfluous.
 C_H2O = sum((8/3*pi*rho_kgm3(1:npre).*r_m(1:npre).^4.*dz_m(1:npre))')'; %calculate C_H2O to the beginning of the depth where we start probing for the beginning of the silicate layer

% Preallocate
C2inds = cell(1);
[C2mean, C2max, C2min, R_sil_mean_m, M_H2O_mean_kg] = deal(zeros(1));
[R_sil_m, rho_sil_kgm3] = deal(zeros(1,nR));
    
    % without a core
if ~Planet.FeCore
    C1 = zeros(1,nR);
    for iz = 1:nR
        C_H2O(iz+1) = C_H2O(iz)+(8/3*pi*rho_kgm3(iz+npre).*r_m(iz+npre).^4.*dz_m(iz+npre));
        R_sil_m(iz) = r_m(iz+npre);
        rho_sil_kgm3(iz) = 3/4/pi*(Planet.M_kg-M_above_kg(iz+npre))./power(R_sil_m(iz),3);
        C1(iz) = C_H2O(iz+1)+8/15*pi*power(R_sil_m(iz),5).*rho_sil_kgm3(iz);
    end
        C2inds = find(C1/MR2>Planet.Cmeasured-Planet.Cuncertainty & C1/MR2<Planet.Cmeasured+Planet.Cuncertainty);
        if isempty(C2inds)
            error(['C/MR2=' num2str(Planet.Cmeasured) ' not found. min:' num2str(min(C1/MR2)) '; max:' num2str(max(C1/MR2))])
        end
        C2mean = round(mean(C2inds));
        C2max = max(C2inds);
        C2min = min(C2inds);
        CMR2mean = C1(C2mean)/MR2;
        R_sil_mean_m = R_sil_m(C2mean);
    R_sil_range_m = R_sil_m(C2min)-R_sil_m(C2max);
        M_H2O_mean_kg = M_above_kg(C2mean);
    R_Fe_mean_m = zeros(1);
   
%     dz_ocean_m_m(~dindsVI & ~dindsV) = Planet.R_m - R_sil_mean_m(~dindsVI &
%     ~dindsV)-zI_m(~dindsVI & ~dindsV); not sure this is right. commenting
%     out on March 8 2018
    rho_Fe_kgm3 = 0;
% =====
else % WITH A CORE
    rho_Fe_kgm3 = Planet.rhoFe*Planet.rhoFeS/(Planet.xFeS*(Planet.rhoFe-Planet.rhoFeS)+Planet.rhoFeS);
    [C2,M_iron_kg,R_Fe_m] = deal(zeros(1,nR));
    for iz = 1:nR
        C_H2O(iz+1) = C_H2O(iz)+abs((8/3*pi*rho_kgm3(iz+npre).*r_m(iz+npre).^4.*dz_m(iz+npre)));
        R_sil_m(iz) = r_m(iz+npre);
        [C2(iz),R_Fe_m(iz)] = CoreSize(Planet.rho_sil_withcore_kgm3,rho_Fe_kgm3,C_H2O(iz+1),M_above_kg(iz+npre),R_sil_m(iz), Planet.M_kg, Planet.R_m);
    end
    
    [R_Fe_mean_m, R_Fe_range_m, R_sil_range_m] = deal(zeros(1));
    r_core_m = zeros(1,Params.nsteps_core);
    C2inds = find(C2/MR2>Planet.Cmeasured-Planet.Cuncertainty & C2/MR2<Planet.Cmeasured+Planet.Cuncertainty);
    if isempty(C2inds)
        error(['C/MR2=' num2str(Planet.Cmeasured) ' not found. min:' num2str(min(C2/MR2)) '; max:' num2str(max(C2/MR2))])
    end
    C2mean = round(mean(C2inds)); % MJS 2021-10-30: This is not a mean value, this operation is taking the mean of indices and then rounding to the nearest index
    C2max = max(C2inds);
    C2min = min(C2inds);
    CMR2mean = C2(C2mean)/MR2;
    R_Fe_mean_m = R_Fe_m(C2mean);
    R_sil_mean_m = R_sil_m(C2mean);
    M_H2O_mean_kg = M_above_kg(C2mean);
    R_Fe_range_m = R_Fe_m(C2max)-R_Fe_m(C2min);
    R_sil_range_m = R_sil_m(C2min)-R_sil_m(C2max);

    r_core_m = linspace(R_Fe_mean_m,0,Params.nsteps_core);

end


%% Print the depths to the various layers
%allocate
zI_m = Planet.zb_outerIce_m;
%zI_m = zb_outerIce_m-Zclath;
[zIII_m, zV_m, zVI_m, indSil] = deal(zeros(1));
indsV = zV_m;
indsVI = zV_m;
indsIII = zV_m;
indsClath = zV_m;
dz_ocean_m_m = Planet.zb_outerIce_m;

% figure out the indices for the tops of the different ice layers
  theseinds = find(phase==1);
if ~isempty(theseinds)
    indsI = theseinds(1);
    zI_m2 = z_m(indsI);
elseif isempty(theseinds) & Params.nsteps_clath>0
    %indsI = theseinds(1);
    zI_m2 = 0;
end

theseinds = find(phase==3);
if ~isempty(theseinds)
    indsIII = theseinds(1);
    zIII_m = z_m(indsIII);
end
theseinds = find(phase==5);
if ~isempty(theseinds)
    indsV = theseinds(1);
    zV_m = z_m(indsV);
end
theseinds = find(phase==6);
if ~isempty(theseinds)
    indsVI = theseinds(1);
    zVI_m = z_m(indsVI);
end

indSil = find(R_sil_m==R_sil_mean_m);

nsteps_mantle = Params.nsteps_mantle;

[dz_ocean_m_m,dzIII_m,dzV_m,dzVI_m,dzClath_m] = deal(zeros(1));
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
if ~isempty(dindsVI) && dindsVI>0
    dzVI_m = Planet.R_m-R_sil_mean_m - zVI_m;
elseif ~isempty(dindsV) && dindsV>0
    dzVI_m = 0;
    dzV_m = RV_m-RVI_m*dindsVI; % this "fix" may introduce an failure condition, but the previous method on the next line was also failing.
%         dzV_m = RV_m-R_sil_mean_m - zVI_m*dindsVI;
elseif ~isempty(dindsIII) && dindsIII>0
    dzV_m = 0;
    dzVI_m = 0;
    dzIII_m = Planet.R_m-R_sil_mean_m - zVI_m;
%   elseif ~isempty(dindsII) && dindsII>0
%       dzV_m = 0;
%       dzVI_m = 0;
%       dzIII_m = Planet.R_m-R_sil_mean_m - zVI_m;
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
[r_Planet_m, sig] = deal(zeros(1,nsteps_tot));
[ocean_thk, ind_Ih, ind_Obot, kmean, ktop] = deal(zeros(1));
Planet.ice_thk = strings(1);
totL = size(r_Planet_m,2);
% Grab the indices for the ices, but only for the portion of the
% interior calculation above the silicate interface
H2Oinds = 1:indSil-1;
indsI = find(phase(H2Oinds)==1);
indsSurf = find(phase(H2Oinds)==1 | phase(H2Oinds)==30);% changes to look for top ice/clath layer
indsClath = find(phase(H2Oinds)==30);
indsLiquid = find(phase(H2Oinds)==0);
indsIII = find(phase(H2Oinds)==3);
indsV = find(phase(H2Oinds)==5);
indsVI = find(phase(H2Oinds)==6);
sig(indsI) = iceSig;
if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
    sig(indsLiquid) = k_S_m(indsLiquid);
else
    sig(indsLiquid) = 0.0;
end
sig(indsIII) = iceSig;
sig(indsV) = iceSig;
sig(indsVI) = iceSig;
sig(indSil:indSil+nsteps_mantle-1) = mantleSig;
if Planet.FeCore
    R_sil_bot = R_Fe_mean_m - (R_Fe_mean_m - R_sil_mean_m)/(nsteps_mantle-1);
    r_mantle_m = linspace(R_sil_mean_m,R_sil_bot,nsteps_mantle);
    sig((nsteps_tot - Params.nsteps_core):end) = coreSig;
    r_Planet_m = [r_m(1:indSil-1) r_mantle_m r_core_m];
else
    r_mantle_m = linspace(R_sil_mean_m,0,nsteps_mantle);
    r_Planet_m = [r_m(1:indSil-1) r_mantle_m];
end
sig(sig==0) = 1e-16; % zeros make the integration unstable
sig(isnan(sig)) = 1e-16;

% Identify lower bound of top ice I/clathrate layer
if ~isempty(indsI)
    ind_Ih = indsI(end);
else
    ind_Ih = indsSurf(end);
end
% Get ocean thickness
if isempty(indsLiquid)
    ocean_thk = 0;
    ind_Obot = ind_Ih;
    ktop = sig(ind_Ih);
    kmean = ktop;
else
    ocean_thk = r_Planet_m(indsLiquid(1)) - r_Planet_m(indsLiquid(end)+1);
    ind_Obot = indsLiquid(end);
    kmean = mean(sig(ind_Ih+1:ind_Obot));
    ktop = sig(ind_Ih+1);
end
disp(['Ocean thickness: ' num2str(ocean_thk*1e-3,'%0.1f km')])


%% Interpolate fewer ocean layers to reduce computational load
if cfg.REDUCED && Params.INCLUDE_ELECTRICAL_CONDUCTIVITY && ocean_thk>0
    nIntL = cfg.nIntL;
    ocStart = totL - indsLiquid(end) + 1;
    ocEnd = totL - indsLiquid(1) + 1;
    nOis = length(indsLiquid);
    interpBds = zeros(1,length(r_Planet_m) - nOis + nIntL);
    interpSig = zeros(1,length(r_Planet_m) - nOis + nIntL);

    r_asc = flip(r_Planet_m,2);
    sig_asc = flip(sig,2);

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
    Planet.boundaries = [interpBds(1:end-1) padBds interpBds(end)];
    Planet.sig = [interpSig(1:end-1) padSig interpSig(end)];
else
    Planet.boundaries = flip(r_Planet_m,2);
    Planet.sig = flip(sig,2);
end

Planet.kmean = kmean(:);
Planet.ktop = ktop(:);

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
%         if CONVECTION_FLAG_I || (isfield(Params,'FORCE_ICEI_CONVECTION') && Params.FORCE_ICEI_CONVECTION == true)
%             %conductive upper layer
%             nConvectIce=n_iceI-nIceIIILithosphere-1;
%             inds_eTBL = find(z_m(1:nConvectIce)<=eTBL_m);
%             Pterm = P_MPa(inds_eTBL)./P_MPa(inds_eTBL(end));
%             T_K(inds_eTBL) = (Tc.^(Pterm)).*(Planet.Tsurf_K.^(1-Pterm));
%             %convective region
%             for iconv = inds_eTBL(end)+1:nConvectIce
%     %             rho = 1000./getVspChoukroun2010(P_MPa(iconv),T_K(iconv-1),2);
%                 rho = getRhoIce(P_MPa(iconv),T_K(iconv-1),1);
%                 try
%                     [Cp,alpha_K] = getCpIce(P_MPa(iconv),T_K(iconv-1),1);
%                 catch
%                     warning('Seafreeze couldn''t get Cp, alpha for ice. Is it installed?');
%                 end
%     %             aK = 1.56e-4;
%                 T_K(iconv) = T_K(iconv-1)+alpha_K*T_K(iconv)/Cp/rho*deltaP*1e6;
%             end
%            % conductive lower layer
%            inds_deltaTBL = find(z_m(1:nConvectIce)>=z_m(nConvectIce)-deltaTBL_m);
%            T_K(inds_deltaTBL) = (Planet.Tb_K.^(P_MPa(inds_deltaTBL)./PbI_MPa)).*(T_K(inds_deltaTBL(1)-1).^(1-P_MPa(inds_deltaTBL)./PbI_MPa));
%            rho_kgm3(inds_deltaTBL) = getRhoIce(P_MPa(inds_deltaTBL),T_K(inds_deltaTBL),1);
%
%            if find(phase>1)
%     %         indVI = find(phase==6);
%     %         Ttop = T_K(indsVI);
%     %         Tbottom = zVI_m
%     %         [Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,CONVECTION_FLAG_I]=...
%     %         ConvectionHPIceKalousova2018(Ttop,,PbI_MPa/2,zb_outerIce_m,g_ms2(1),2);
%            end
%
%         %else % We take care of this case well above in order to be able
%               % to exit sooner if cfg.SKIP_PROFILES=1
%         end
%     else
%         Q_Wm2 = Planet.Qmantle_Wm2(1);
%     end
    
%end
%% Calculate Ice velocities
if Params.CALC_NEW_SOUNDSPEEDS
    disp('Computing sound speeds')
%     velsIce = iceVelsGagnon1990(P_MPa,T_K);
    [vfluid_kms,Ksfluid_GPa] = deal(zeros(1,nsteps));
    tic
    out=Helgerud_sI(P_MPa', T_K');
    velsIce.Vclathl_kms = 1e-3*out.Vp;
    velsIce.Vclatht_kms = 1e-3*out.Vs;
    velsIce.Ksclath_GPa = 1e-3*out.K;
    velsIce.Gsclath_GPa = 1e-3*out.shear;

    out=SeaFreeze([P_MPa' T_K'],'Ih');
    velsIce.VIl_kms = 1e-3*out.Vp;
    velsIce.VIt_kms = 1e-3*out.Vs;
    velsIce.KsI_GPa = 1e-3*out.Ks;
    velsIce.GsI_GPa = 1e-3*out.shear;

    out=SeaFreeze([P_MPa' T_K'],'II');
    velsIce.VIIl_kms = 1e-3*out.Vp;
    velsIce.VIIt_kms = 1e-3*out.Vs;
    velsIce.KsII_GPa = 1e-3*out.Ks;
    velsIce.GsII_GPa = 1e-3*out.shear;

    out=SeaFreeze([P_MPa' T_K'],'III');
    velsIce.VIIIl_kms = 1e-3*out.Vp;
    velsIce.VIIIt_kms = 1e-3*out.Vs;
    velsIce.KsIII_GPa = 1e-3*out.Ks;
    velsIce.GsIII_GPa = 1e-3*out.shear;

    out=SeaFreeze([P_MPa' T_K'],'V');
    velsIce.VVl_kms = 1e-3*out.Vp;
    velsIce.VVt_kms = 1e-3*out.Vs;
    velsIce.KsV_GPa = 1e-3*out.Ks;
    velsIce.GsV_GPa = 1e-3*out.shear;

    out=SeaFreeze([P_MPa' T_K'],'VI');
    velsIce.VVIl_kms = 1e-3*out.Vp;
    velsIce.VVIt_kms = 1e-3*out.Vs;
    velsIce.KsVI_GPa = 1e-3*out.Ks;
    velsIce.GsVI_GPa = 1e-3*out.shear;


    if POROUS_ICE
        % only affects top ice/clath layers
        % correction for porosity
        por_in.p=P_MPa(1:(n_clath+n_iceI))*1e-3;
        por_in.t = T_K(1:(n_clath+n_iceI));
        por_in.den = rho_kgm3(1:(n_clath+n_iceI));
        try
            por_in.phi_surface=Planet.phi_surface;
        end

        por_in.vp = velsIce.Vclathl_kms(1:(n_clath+n_iceI));
        por_in.vs = velsIce.Vclatht_kms(1:(n_clath+n_iceI));
        por_vel_out = porosity_correction_ice(por_in,por);
        velsIce.Vclathl_kms(1:(n_clath+n_iceI)) = por_vel_out.vp;
        velsIce.Vclatht_kms(1:(n_clath+n_iceI)) = por_vel_out.vs;


        por_in.vp = velsIce.VIl_kms(1:(n_clath+n_iceI));
        por_in.vs = velsIce.VIt_kms(1:(n_clath+n_iceI));
        por_vel_out = porosity_correction_ice(por_in,por);
        velsIce.VIl_kms(1:(n_clath+n_iceI)) = por_vel_out.vp;
        velsIce.VIt_kms(1:(n_clath+n_iceI)) = por_vel_out.vs;
    end

    ir = find(phase==0);
    vfluid_kms(ir) = fluidSoundSpeeds(P_MPa(ir),T_K(ir),wo,Planet.Ocean.comp);
    ireal = find(~isnan(vfluid_kms));
    if length(ireal)<length(ir)
        if ~strcmp(Planet.Ocean.comp,'Seawater')
            disp('WARNING: extrapolating fluid sound speeds, but this is only valid for Seawater above 120 MPa!')
        end
        vfluid_kms(ir) = interp1(P_MPa(ireal),vfluid_kms(ireal),P_MPa(ir),'linear','extrap');
    end
%         for ir = find(phase==0)
%             vfluid_kms(ir) = fluidSoundSpeeds(P_MPa(ir),T_K(ir),wo,Planet.Ocean.comp);
%         end
    toc
    Ksfluid_GPa(ir)=1e-3*rho_kgm3(ir).*vfluid_kms(ir).^2;
    vfluid_kms(vfluid_kms==0)=NaN;
    Ksfluid_GPa(vfluid_kms==0)=NaN;
    velsIce.VIl_kms(phase~=1) =NaN;
    velsIce.VIt_kms(phase~=1) =NaN;
    velsIce.VIIl_kms(phase~=2) =NaN;
    velsIce.VIIt_kms(phase~=2) =NaN;
    velsIce.VIIIl_kms(phase~=3) =NaN;
    velsIce.VIIIt_kms(phase~=3) =NaN;
    velsIce.VVl_kms(phase~=5) =NaN;
    velsIce.VVt_kms(phase~=5) =NaN;
    velsIce.VVIl_kms(phase~=6) =NaN;
    velsIce.VVIt_kms(phase~=6) =NaN;
    velsIce.Vclathl_kms(phase~=30) =NaN;
    velsIce.Vclatht_kms(phase~=30) =NaN;
    save(fullfile([datpath savefile '_Vels' '_pp']),'Ksfluid_GPa','velsIce','vfluid_kms');
else
    try	
        load(fullfile([datpath savefile '_Vels' '_pp']),'Ksfluid_GPa','velsIce','vfluid_kms');	
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
    

mean_rho_ocean = mean(rho_kgm3(n_iceI+1:indSil));
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

[Pcore_MPa,Tcore_K,rhocore_kgm3,M_above_core,VP_core_kms,VS_core_kms,g_ms2_core,Ks_iron_GPa,G_iron_GPa] = deal(zeros(1,Params.nsteps_core));
clear g_ms2_core M_above_core
    
[Pmantle_MPa,rho_mantle_kgm3,M_above_mantle,VP_mantle_kms,VS_mantle_kms,g_ms2_sil] = deal(zeros(1,nsteps_mantle));
Tmantle_K = T_K(indSil);
r_mantle_m = linspace(R_sil_mean_m,R_Fe_mean_m,nsteps_mantle);

g_ms2_sil = g_ms2(C2mean)*ones(1,nsteps_mantle);
MantleHeat = Planet.Qmantle_Wm2*4*pi*Planet.R_m^2+Planet.QHmantle;
rhom_rough = 3000;
alpha_rough = 0.2e-4;
Cp_rough = 2e6;
Kappa_rough = 1e-6;
nu_mantle = 1e21; % mantle viscosity in Pa S, a common number for Earth?s mantle
DeltaT = 800; % temperature differential in K

% cold case
if Planet.FeCore
    Tmantle_K = conductiveMantleTemperature(r_mantle_m,R_Fe_mean_m,R_sil_mean_m,Planet.kr_mantle,Planet.rho_sil_withcore_kgm3,Tmantle_K,MantleHeat,0);
    MantleMass = 4/3*pi*(R_sil_mean_m^3-R_Fe_mean_m^3);
    Dm = R_sil_mean_m-R_Fe_mean_m;
else
    Tmantle_K = conductiveMantleTemperature(r_mantle_m,0,R_sil_mean_m,Planet.kr_mantle,rho_sil_kgm3(C2mean),Tmantle_K,MantleHeat,0);
    MantleMass = 4/3*pi*(R_sil_mean_m^3);
    Dm = R_sil_mean_m;
end    

Ra_mantle = g_ms2_sil(1)*rhom_rough*alpha_rough*MantleHeat/MantleMass*Dm^5/nu_mantle/Kappa_rough^2/Cp_rough; % bunge 1997
% parmentier and sotin 1994
Ra_mantle = alpha_rough*rhom_rough*g_ms2_sil(1)*DeltaT*Dm^3/Kappa_rough/nu_mantle;% as per Hussmann and Spohn 2004

% COMMON FUNCTIONALITY HERE IS TO PROPAGATE PRESSURE, GRAVITY AND
% DENSITY DOWNWARD FOR A CONDUCTIVELY COOLING LAYER WITH NO SOLID STATE
% CONVECTION
%     g_ms2_sil(1) = g_ms2(C2mean);
M_above_mantle(1) = M_above_kg(C2mean);
Pmantle_MPa(1) = P_MPa(C2mean);
rhofn_kgm3 = mantle.rho_fn;
VPfn_kms = mantle.vp_fn;
VSfn_kms = mantle.vs_fn;
rho_mantle_kgm3(1) = rhofn_kgm3(Pmantle_MPa(1),Tmantle_K(1));
VP_mantle_kms(1) = VPfn_kms(Pmantle_MPa(1),Tmantle_K(1));
VS_mantle_kms(1) = VSfn_kms(Pmantle_MPa(1),Tmantle_K(1));
Cpfn_mantle_JkgK = mantle.cp_fn;
alphafn_mantle_pK = mantle.cp_fn;

for ij = 2:nsteps_mantle
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
    por_in.p=Pmantle_MPa*1e-3;
    por_in.t = Tmantle_K;
    por_in.den = rho_mantle_kgm3;
    por_in.vp = VP_mantle_kms;
    por_in.vs = VS_mantle_kms;
    if isfield(Planet,'phi_surface')
        por_out = porosity_correction(por_in,Planet.phi_surface);
    else
        por_out = porosity_correction(por_in);
    end
    permeability = por_out.per;
    rho_mantle_kgm3 = por_out.den;
    VP_mantle_kms = por_out.vp;
    VS_mantle_kms = por_out.vs;
    phi_mantle_frac = por_out.por;
    disp(['Average Porosity: ' num2str(mean(phi_mantle_frac))])
    disp(['Porosity: ' num2str(phi_mantle_frac)])
    %recalculate mass above in light of reduced density. This should
    %really be done in a recursive way that acknowledges the reduced
    %overburden pressure.
    for ij = 2:nsteps_mantle
%         if r_mantle_m(ij-1)>0.3*Planet.R_m
%             g_ms2_sil(ij)=(Gg*(Planet.M_kg-M_above_mantle(ij-1))/r_mantle_m(ij-1)^2);
%         else
%             g_ms2_sil(ij) = g_ms2_sil(ij-1);
%         end
%         Pmantle_MPa(ij) = Pmantle_MPa(ij-1)+1e-6*(rho_mantle_kgm3(ij-1))*g_ms2_sil(ij)*(r_mantle_m(ij-1)-r_mantle_m(ij));
      M_above_mantle(ij) = M_above_mantle(ij-1)+4/3*pi*(r_mantle_m(ij-1)^3-r_mantle_m(ij)^3)*rho_mantle_kgm3(ij-1);
    end
end

Cp_mantle_JkgK = Cpfn_mantle_JkgK(Pmantle_MPa,Tmantle_K);
alpha_mantle_pK = alphafn_mantle_pK(Pmantle_MPa,Tmantle_K);

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
%     g_ms2_core(1) = g_ms2_sil(nsteps_mantle); % depth dependent
%     gravity was used previously, but this can result in negative gravity
%     approaching the center of the core due to imperfect matching of the
%     moment of inertia.
%     M_above_core(1) = M_above_mantle(nsteps_mantle);
if Planet.FeCore
    M_above_core(1) = M_above_mantle(nsteps_mantle);
    g_ms2_core=(Gg*(Planet.M_kg-M_above_core)/r_core_m(1)^2)*ones(1,Params.nsteps_core);

    Tcore_K = linspace(Tmantle_K(end),1.01*Tmantle_K(end),Params.nsteps_core);
    %Pcore_MPa(1) = 2*Pmantle_MPa(nsteps_mantle) - Pmantle_MPa(nsteps_mantle-1);
    Pcore_MPa(1) = Pmantle_MPa(nsteps_mantle);
    if isfield(Seismic,'coreEOS')
        core = loadMantleEOS(Seismic.coreEOS);
        rhocore_kgm3(1) = core.rho_fn(Pcore_MPa(1),Tcore_K(1));

        for ij = 2:Params.nsteps_core
    %          g_ms2_core(ij)=(Gg*(Planet.M_kg-M_above_core(ij-1))/r_core_m(ij-1)^2);
    %          Pcore_MPa(ij) = Pcore_MPa(ij-1)+1e-6*rhocore_kgm3(ij-1)*g_ms2_core(ij)*(r_core_m(ij-1)-r_core_m(ij));
             Pcore_MPa(ij) = Pcore_MPa(ij-1)+1e-6*rhocore_kgm3(ij-1)*g_ms2_core(ij-1)*(r_core_m(ij-1)-r_core_m(ij));
             rhocore_kgm3(ij) = core.rho_fn(Pcore_MPa(ij),Tcore_K(ij));
%              rhocore_kgm3(ij) = rhocore_kgm3(ij-1)*(1+1./(1e-3*Pcore_MPa(ij)*Ks_iron_GPa(ij)));
    %          M_above_core(ij) = M_above_core(ij-1)+4/3*pi*(r_core_m(ij-1)^3-r_core_m(ij)^3)*rhocore_kgm3(ij-1);
        end
        Ks_iron_GPa = bar2GPa*core.Ks_fn(Pcore_MPa,Tcore_K);
        G_iron_GPa = bar2GPa*core.Gs_fn(Pcore_MPa,Tcore_K);
        VS_core_kms = core.vs_fn(Pcore_MPa,Tcore_K);
        VP_core_kms = core.vp_fn(Pcore_MPa,Tcore_K);

        Cpfn_core_JkgK = core.cp_fn;
        alphafn_core_pK = core.alpha_fn;
        Cp_core_JkgK = Cpfn_core_JkgK(Pcore_MPa,Tcore_K);
        alpha_core_pK = alphafn_core_pK(Pcore_MPa,Tcore_K);
    else
        Ks_iron_GPa(1) = Kso_iron_GPa+1e-3*Pcore_MPa(1)*dKsdP_iron+1e-9*Tcore_K(1)*dKsdT_PaK;
        G_iron_GPa(1) = Go_iron_GPa+1e-3*Pcore_MPa(1)*dGdP_iron+1e-9*Tcore_K(1)*dGdT_iron_PaK;
        rhocore_kgm3 = rhoo_iron_kgm3;

        for ij = 2:Params.nsteps_core
    %          g_ms2_core(ij)=(Gg*(Planet.M_kg-M_above_core(ij-1))/r_core_m(ij-1)^2);
    %          Pcore_MPa(ij) = Pcore_MPa(ij-1)+1e-6*rhocore_kgm3(ij-1)*g_ms2_core(ij)*(r_core_m(ij-1)-r_core_m(ij));
             Pcore_MPa(ij) = Pcore_MPa(ij-1)+1e-6*rhocore_kgm3(ij-1)*g_ms2_core*(r_core_m(ij-1)-r_core_m(ij));
             Ks_iron_GPa(ij) = Kso_iron_GPa+1e-3*Pcore_MPa(ij)*dKsdP_iron+1e-9*Tcore_K(ij)*dKsdT_PaK;
             G_iron_GPa(ij) = Go_iron_GPa+1e-3*Pcore_MPa(ij)*dGdP_iron+1e-9*Tcore_K(ij)*dGdT_iron_PaK;
%              rhocore_kgm3(ij) = rhocore_kgm3(ij-1)*(1+1./(1e-3*Pcore_MPa(ij)*Ks_iron_GPa(ij)));
    %          M_above_core(ij) = M_above_core(ij-1)+4/3*pi*(r_core_m(ij-1)^3-r_core_m(ij)^3)*rhocore_kgm3(ij-1);
        end
        VS_core_kms = 1e-3*sqrt(G_iron_GPa*1e9./rhocore_kgm3);
        VP_core_kms = 1e-3*sqrt(Ks_iron_GPa*1e9./rhocore_kgm3+4/3*(VS_core_kms*1e3).^2);
    end

end

%     kluge to fix NaN's in the mantle densities and sound speeds
rho_mantle_kgm3(isnan(rho_mantle_kgm3))=0;
VS_mantle_kms(isnan(VS_mantle_kms))=0;
VP_mantle_kms(isnan(VP_mantle_kms))=0;


interior.g_ms2 = g_ms2_sil;
interior.r_mantle_m = r_mantle_m;
interior.rho_mantle_kgm3 = rho_mantle_kgm3;
interior.Tmantle_K = Tmantle_K;
interior.Pmantle_MPa = Pmantle_MPa;
interior.VS_mantle_kms = interp1nan(Pmantle_MPa,VS_mantle_kms);
interior.VP_mantle_kms = interp1nan(Pmantle_MPa,VP_mantle_kms);
if POROUS
    interior.permeability = permeability;
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
            interior.(['fluid_' mfnames{im}]) = mfluids.(mfnames{im+1})(Pmantle_MPa,Tmantle_K);
        end
    end
    for im = 1:length(mpnames) % count only the names, not the interpolating functions or p or t
        if ~endsWith(mpnames{im},'_fn') && ~(strcmp(mpnames{im},'p') || strcmp(mpnames{im},'t'))
            interior.(['wt_' mpnames{im}]) = mphases.(mpnames{im+1})(Pmantle_MPa,Tmantle_K); %wt pct of each constituent
        end
    end
    masscheck = 0;
    for im = 1:length(mvnames) % count only the names, not the interpolating functions or p or t
        if ~endsWith(mvnames{im},'_fn') && ~(strcmp(mvnames{im},'p') || strcmp(mvnames{im},'t'))
            disp(mvnames{im})
            interior.(['vol_' mvnames{im}]) = mvolumes.(mvnames{im+1})(Pmantle_MPa,Tmantle_K); %vol pct of each constituent
            interior.(['rho_' mvnames{im}]) = interior.(['wt_' mpnames{im}])./interior.(['vol_' mvnames{im}]).*rho_mantle_kgm3; %dens fraction of each constituent
            interior.(['mass_' mvnames{im}]) = 4/3*pi*interior.(['rho_' mvnames{im}]).*-gradient(r_mantle_m.^3); %mass in each spherical layer of the mantle
        end
    end
end
if isfield(mantle,'Ks_fn')
interior.Ks_GPa = bar2GPa*mantle.Ks_fn(Pmantle_MPa,Tmantle_K);
interior.Ks_GPa = interp1nan(Pmantle_MPa,interior.Ks_GPa);
interior.Gs_GPa = bar2GPa*mantle.Gs_fn(Pmantle_MPa,Tmantle_K);
interior.Gs_GPa = interp1nan(Pmantle_MPa,interior.Gs_GPa);
end
    % compute seismic attenuation, as per C2006,
% Rb = 8.314462175; % J/K/mol
Tm_K = 273.15+Tm_p_Hirschmann2000(1e-3*Pmantle_MPa); %input in GPa
Hp_mantle = Seismic.g_aniso_mantle*Tm_K;
interior.QS_overfgamma = Seismic.B_aniso_mantle*exp(Seismic.gamma_aniso_mantle*Hp_mantle./Tmantle_K);

%% Construct output for seismic modeling:
    % ice                           ocean                                            mantle         core
% [g_Planet_ms2,P_Planet_MPa,T_Planet_K,r_Planet_m,rho_pPlanet_kgm3,VP_Planet_kms,VS_Planet_kms,QS_overfgamma_Planet,k_S_mPlanet]=...
%     deal(zeros(nTbs,indSil(1)-1+nsteps_mantle(1)+Params.nsteps_core));
%% Plot settings
ymax = 1.05*Planet.R_m*1e-3;

if Planet.FeCore
[VP_Planet_kms,VS_Planet_kms,Ks_Planet_GPa,Gs_Planet_GPa,QS_overfgamma_Planet,k_S_mPlanet,phasePlanet] = ...
    deal(nan(1,length([g_ms2(1:indSil(1)-1) interior.g_ms2 g_ms2_core(1,:)])));
else
[VP_Planet_kms,VS_Planet_kms,Ks_Planet_GPa,Gs_Planet_GPa,QS_overfgamma_Planet,k_S_mPlanet,phasePlanet] = ...
    deal(nan(1,length([g_ms2(1:indSil(1)-1) interior.g_ms2])));
end

% Grab the indices for the ices, but only for the portion of the
% interior calculation above the silicate interface
H2Oinds = 1:indSil-1;
indsI = find(phase(H2Oinds)==1);
indsLiquid = find(phase(H2Oinds)==0);
indsIII = find(phase(H2Oinds)==3);
indsV = find(phase(H2Oinds)==5);
indsVI = find(phase(H2Oinds)==6);
indsClath = find(phase(H2Oinds)==30);

phasePlanet(indsI) = 1;
phasePlanet(indsLiquid) = 0;
phasePlanet(indsIII) = 3;
phasePlanet(indsV) = 5;
phasePlanet(indsVI) = 6;
phasePlanet(indsClath) = 30;

phasePlanet(indSil:indSil+nsteps_mantle-1)=50; %reserve 50-99 for different mantle types

if ~Planet.FeCore
    [thisgcore,thisPcore,thisTcore,thisrcore,thisrhocore,thisVPcore,thisVScore,thisQScore,thisKScore,thisGScore,thiskScore]=deal([]); %placeholders to avoid errors
else
    thisgcore = g_ms2_core;
    thisPcore = Pcore_MPa;
    thisTcore = Tcore_K;
    thisrcore = r_core_m;
    thisrhocore = rhocore_kgm3;
    thisVPcore = VP_core_kms;
    thisVScore = VS_core_kms;
    thisQScore = Seismic.QScore*ones(1,Params.nsteps_core);
    thiskScore = 0*Seismic.QScore*ones(1,Params.nsteps_core);
    thisKScore = Ks_iron_GPa;
    thisGScore = G_iron_GPa;
    phasePlanet(end-Params.nsteps_core+1:end) = 100;
end
g_Planet_ms2 = [g_ms2(1:indSil-1) interior.g_ms2 thisgcore];
P_Planet_MPa = [P_MPa(1:indSil-1) interior.Pmantle_MPa thisPcore];
T_Planet_K = [T_K(1:indSil-1) interior.Tmantle_K thisTcore];
r_Planet_m = [r_m(1:indSil-1)    interior.r_mantle_m thisrcore];
rho_pPlanet_kgm3 = [rho_kgm3(1:indSil-1) interior.rho_mantle_kgm3 thisrhocore];
ocean_thk = r_Planet_m(indsLiquid(1)) - r_Planet_m(indsLiquid(end)+1);
ind_Obot = indsLiquid(end);

VP_Planet_kms(indsI) = velsIce.VIl_kms(indsI);
VP_Planet_kms(indsLiquid) = vfluid_kms(indsLiquid);
VP_Planet_kms(indsIII) = velsIce.VIIIl_kms(indsIII);
VP_Planet_kms(indsV) = velsIce.VVl_kms(indsV);
VP_Planet_kms(indsClath) = velsIce.Vclathl_kms(indsClath);
VP_Planet_kms(indsVI) = velsIce.VVIl_kms(indsVI);
VP_Planet_kms(indSil:indSil+length(interior.VP_mantle_kms)-1) = interior.VP_mantle_kms;
if Planet.FeCore
    start_core = indSil+length(interior.VP_mantle_kms);
    VP_Planet_kms(start_core:start_core-1+length(thisVPcore)) = thisVPcore;
end

%     VS_Planet_kms = [velsIce.VIt_kms(1:n_iceI) 0*vfluid_kms(indsLiquid) ...
%         velsIce.VIIIt_kms(indsIII) velsIce.VVt_kms(indsV) velsIce.VVIt_kms(indsVI) ...
%         interior.VS_mantle_kms thisVScore];
VS_Planet_kms(indsI) = velsIce.VIt_kms(indsI);
VS_Planet_kms(indsLiquid) = 0*vfluid_kms(indsLiquid);
VS_Planet_kms(indsIII) = velsIce.VIIIt_kms(indsIII);
VS_Planet_kms(indsV) = velsIce.VVt_kms(indsV);
VS_Planet_kms(indsClath) = velsIce.Vclatht_kms(indsClath);
VS_Planet_kms(indsVI) = velsIce.VVIt_kms(indsVI);
VS_Planet_kms(indSil:indSil+length(interior.VS_mantle_kms)-1) = interior.VS_mantle_kms;
if Planet.FeCore
    VS_Planet_kms(start_core:start_core-1+length(thisVPcore)) = thisVScore;
end

if isfield(mantle,'Ks_fn')
%     Ks_Planet_GPa=[velsIce.KsI_GPa(1:n_iceI) Ksfluid_GPa(indsLiquid) ...
%         velsIce.KsIII_GPa(indsIII) velsIce.KsV_GPa(indsV) velsIce.KsVI_GPa(indsVI) ...
%         interior.Ks_GPa thisKScore];
    Ks_Planet_GPa(indsI) = velsIce.KsI_GPa(indsI);
    Ks_Planet_GPa(indsLiquid) = Ksfluid_GPa(indsLiquid);
    Ks_Planet_GPa(indsIII) = velsIce.KsIII_GPa(indsIII);
    Ks_Planet_GPa(indsV) = velsIce.KsV_GPa(indsV);
    Ks_Planet_GPa(indsClath) = velsIce.Ksclath_GPa(indsClath);
    Ks_Planet_GPa(indsVI) = velsIce.KsVI_GPa(indsVI);
    Ks_Planet_GPa(indSil:indSil+length(interior.VS_mantle_kms)-1) = interior.Ks_GPa;
    if Planet.FeCore
        Ks_Planet_GPa(start_core:start_core-1+length(thisVPcore)) = thisKScore;
    end

%     Gs_Planet_GPa=[velsIce.GsI_GPa(1:n_iceI) 0*Ksfluid_GPa(indsLiquid) ...
%         velsIce.GsIII_GPa(indsIII) velsIce.GsV_GPa(indsV) velsIce.GsVI_GPa(indsVI) ...
%         interior.Gs_GPa thisGScore];
    Gs_Planet_GPa(indsI) = velsIce.GsI_GPa(indsI);
    Gs_Planet_GPa(indsLiquid) = 0*Ksfluid_GPa(indsLiquid);
    Gs_Planet_GPa(indsIII) = velsIce.GsIII_GPa(indsIII);
    Gs_Planet_GPa(indsV) = velsIce.GsV_GPa(indsV);
    Gs_Planet_GPa(indsClath) = velsIce.Gsclath_GPa(indsClath);
    Gs_Planet_GPa(indsVI) = velsIce.GsVI_GPa(indsVI);
    Gs_Planet_GPa(indSil:indSil+length(interior.VS_mantle_kms)-1) = interior.Gs_GPa;
    if Planet.FeCore
        Gs_Planet_GPa(start_core:start_core-1+length(thisVPcore)) = thisGScore;
    end

end
%     QS_overfgamma_Planet = [QS_overfgamma_iceI 0*vfluid_kms(indsLiquid)...
%         QS_overfgamma_iceIII(indsIII) QS_overfgamma_iceV(indsV) QS_overfgamma_iceVI(indsVI) ...
%         interior.QS_overfgamma thisQScore];

%% attenuation in ice
Hp_iceI = Seismic.g_aniso_iceI*Planet.Tb_K;
QS_overfgamma_iceI = Seismic.B_aniso_iceI*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceI./T_K(indsI))/Seismic.LOW_ICE_Q;
try
    Tthis=T_K(indsClath);
    Hp_clath = Seismic.g_aniso_clath*Planet.Tb_K;
    QS_overfgamma_clath= Seismic.B_aniso_clath*....
        exp(Seismic.gamma_aniso_clath*Hp_clath./Tthis)/Seismic.LOW_ICE_Q;
catch
    Hp_clath = Seismic.g_aniso_iceI*Planet.Tb_K;
    QS_overfgamma_clath = Seismic.B_aniso_iceI*....
        exp(Seismic.gamma_aniso_iceI*Hp_clath./T_K(indsClath))/Seismic.LOW_ICE_Q;
end
Tthis=T_K(indsIII);
Hp_iceIII = Seismic.g_aniso_iceIII*max(T_K(indsIII));
QS_overfgamma_iceIII = Seismic.B_aniso_iceIII*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceIII./Tthis)/Seismic.LOW_ICE_Q;
Tthis=T_K(indsV);
Hp_iceV = Seismic.g_aniso_iceV*max(T_K(indsV));
QS_overfgamma_iceV = Seismic.B_aniso_iceV*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceV./Tthis)/Seismic.LOW_ICE_Q;
Tthis=T_K(indsVI);
Hp_iceVI = Seismic.g_aniso_iceVI*max(T_K(indsVI));
QS_overfgamma_iceVI = Seismic.B_aniso_iceVI*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceVI./Tthis)/Seismic.LOW_ICE_Q;

QS_overfgamma_Planet(indsI) = QS_overfgamma_iceI;
QS_overfgamma_Planet(indsLiquid) = 0*vfluid_kms(indsLiquid);
QS_overfgamma_Planet(indsIII) = QS_overfgamma_iceIII;
QS_overfgamma_Planet(indsV) = QS_overfgamma_iceV;
QS_overfgamma_Planet(indsVI) = QS_overfgamma_iceVI;
QS_overfgamma_Planet(indsClath) = QS_overfgamma_clath;
QS_overfgamma_Planet(indSil:indSil+length(interior.VS_mantle_kms)-1) = interior.QS_overfgamma;
if Planet.FeCore
    QS_overfgamma_Planet(start_core:start_core-1+length(thisVPcore)) = thisQScore;
end

%%save the data to a text file
Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [P_Planet_MPa', T_Planet_K', r_Planet_m'*1e-3, rho_pPlanet_kgm3',VP_Planet_kms',VS_Planet_kms',QS_overfgamma_Planet' Ks_Planet_GPa' Gs_Planet_GPa' g_Planet_ms2' phasePlanet'];
header = sprintf('%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',...
    'P (MPa)','T (K)','r (km)','rho (kg m-3)','VP (km s-1)','VS (km s-1)','QS/gamma','KS (GPa)','GS (GPa)','g (m/s2)','phase');
if cfg.CONDUCT %11/8/19 added g and phase. This changes the column number of electrical conductivity
    k_S_mPlanet(indsI) = 0*QS_overfgamma_iceI;
    k_S_mPlanet(indsLiquid) = k_S_m(indsLiquid);
    k_S_mPlanet(indsIII) = 0*QS_overfgamma_iceIII;
    k_S_mPlanet(indsV) = 0*QS_overfgamma_iceV;
    k_S_mPlanet(indsClath) = 0*QS_overfgamma_clath;
    k_S_mPlanet(indsVI) = 0*QS_overfgamma_iceVI;
    k_S_mPlanet(indSil:indSil+length(interior.VS_mantle_kms)-1) = 0*interior.QS_overfgamma;
    if Planet.FeCore
        k_S_mPlanet(start_core:start_core-1+length(thisVPcore)) = 0*thisQScore;
    end
    Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma k_S_mPlanet' ];
    header = sprintf('%s%s\t%s\t',header,'k (S m-1)');

    kmean = mean(k_S_mPlanet(indsLiquid));
    ktop = k_S_mPlanet(indsLiquid(1));
end
if isfield(Planet,'Clathrate'); clathID = ['_Zclath' num2str(Zclath./1000,'%2.0f') 'km']; else; clathID = ''; end
thissavestr = [savefile '_Zb' strLow num2str(Planet.zb_outerIce_m./1000,'%2.0f') 'km' clathID ];
dlmwrite(fullfile([datpath thissavestr '_pp' '.txt']),header,'delimiter','');
dlmwrite(fullfile([datpath thissavestr '_pp' '.txt']),Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma,...
    'delimiter','\t',...
    'precision','%3.5e',...
    '-append');
dlmwrite(fullfile([datpath thissavestr '_mantle_pp.dat']),[10*interior.Pmantle_MPa' interior.Tmantle_K'],'delimiter','\t');

if isfield(Params,'TauP')
    taup_model=create_Taup([datpath thissavestr '.txt']);
end

%%
% MJS 2021-10-31: This block has several loops over iT inside of loops
% over iT. There's a lot that could be done to streamline it.
rhommodel = mean(interior.rho_mantle_kgm3);
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

    dV = -gradient([r_m(1:indSil(iT2)) interior.r_mantle_m].^3);
    Planet.Mcomputed_kg = 4/3*pi*sum(dV.*[rho_kgm3(1:indSil) interior.rho_mantle_kgm3]);
    % Preallocate

    disp(['Planet Mass (kg): ' num2str(Planet.M_kg)])
    disp(['Computed Mass (kg): ' num2str(Planet.Mcomputed_kg,' %0.4g')]);
    disp(['Input C/MR2: ' num2str(Planet.Cmeasured)])
    if Planet.FeCore
        disp(['Computed C/MR2 for Tb=' num2str(Planet.Tb_K) 'K: ' num2str(C2(C2mean)/MR2) ' (neighboring values: ' num2str(C2(C2mean+1)/MR2) '; ' num2str(C2(C2mean-1)/MR2) ')'])
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
        disp(['Computed C/MR2 for Tb=' num2str(Planet.Tb_K) 'K: ' num2str(C1(C2mean)/MR2) '  (neighboring values: ' num2str(C1(C2mean+1)/MR2) '; ' num2str(C1(C2mean-1)/MR2) ')'])
        rhom = mean(rho_sil_kgm3(C2mean));
        rhorockstr = getTableStr(Planet.Tb_K,round(rhom));
        rhorockmstr = getTableStr(Planet.Tb_K,round(rhommodel));
        if exist('concstr','var')
            disp([concstr{1} '&$\rho_{rock}$ (kg m$^{-3}$)' rhorockstr '\\']);
            disp([concstr{2} '&$\rho_{rock,model}$ (kg m$^{-3}$)' rhorockmstr '\\']);
        else
            disp(['&$\rho_{rock}$ (kg m$^{-3}$)&' rhorockstr '\\']);
            disp(['&$\rho_{rock,model}$ (kg m$^{-3}$)' rhorockmstr '\\']);
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
        meanphivec = 100*mean(phi_mantle_frac);
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

    savestr = [savebase Planet.Ocean.comp ...
        '_' num2str(round(wtPpt)) 'WtPpt' minEOS porIceStr ...
        '_Tb' num2str(round(Planet.Tb_K,3),'%.3f') 'K' ];

    if Planet.FeCore
        R_Fe_trade_m = R_Fe_m(C2inds);
        rho_sil_trade_kgm3 = zeros(size(R_sil_m(C2inds)));
        rho_sil_mean_kgm3 = Planet.rho_sil_withcore_kgm3;
        phi_core_frac = zeros(1, Params.nsteps_core);
        rho_core_mean_kgm3 = rho_Fe_kgm3;
    else
        R_Fe_trade_m = zeros(size(R_sil_m(C2inds)));
        rho_sil_trade_kgm3 = rho_sil_kgm3(C2inds);
        rho_sil_mean_kgm3 = rhommodel;
        phi_core_frac = zeros(1, 0);
        rho_core_mean_kgm3 = 0;
        R_Fe_range_m = 0;
        Cp_core_JkgK = zeros(1, 0);
        alpha_core_pK = zeros(1, 0);
    end
    R_sil_trade_m = R_sil_m(C2inds);

    if ~Planet.POROUS_ROCK
        phi_mantle_frac = zeros(size(r_mantle_m));
        permeability = zeros(1,Params.nsteps_mantle,5);
    end

    nsteps_hydro = indSil - 1;
    phi_hydro_frac = phi_hydro_frac(1:nsteps_hydro);
    Cp_Planet_JkgK = [Cp(1:nsteps_hydro) Cp_mantle_JkgK Cp_core_JkgK];
    alpha_Planet_pK = [alpha_K(1:nsteps_hydro) alpha_mantle_pK alpha_core_pK];
    phi_Planet_frac = [phi_hydro_frac phi_mantle_frac phi_core_frac];

    writeToDisk([datpath savestr], Planet.Ocean.comp, Planet.FeCore, Planet.Ocean.w_ocean_pct, Planet.Tb_K, ...
        Planet.zb_outerIce_m/1e3, Planet.zClath_m, Pb_MPa, PbI_MPa, deltaP, CMR2mean, ...
        Planet.Qmantle_Wm2, Planet.phi_surface, R_sil_mean_m, R_sil_range_m, rho_sil_mean_kgm3, R_Fe_mean_m, R_Fe_range_m, ...
        rho_core_mean_kgm3, Params.nsteps_clath, Params.nsteps_iceI, nIceIIILithosphere, nIceVLithosphere, ...
        nsteps_hydro, Params.nsteps_mantle, Params.nsteps_core, ...
        R_sil_trade_m, R_Fe_trade_m, rho_sil_trade_kgm3, permeability, ...
        P_Planet_MPa, T_Planet_K, r_Planet_m, phasePlanet, ...
        rho_pPlanet_kgm3, Cp_Planet_JkgK, alpha_Planet_pK, g_Planet_ms2, ...
        phi_Planet_frac, k_S_mPlanet, VP_Planet_kms, VS_Planet_kms, QS_overfgamma_Planet, ...
        Ks_Planet_GPa, Gs_Planet_GPa);
    %save(fullfile([datpath savefile '_pp' num2str]),'P_MPa','Pb_MPa','PbI_MPa','nIceIIILitho','T_K','Tb_K','phase','deltaP','wo','rho_kgm3','Cp','nsteps','n_clath','n_iceI','n_ocean','max_clath_depth'); % save the progress at each step
        

    else
        npre = 0;

        savestr = [savebase Planet.Ocean.comp ...
            '_' num2str(round(wtPpt)) 'WtPpt' minEOS porIceStr ...
            '_Tb' num2str(round(Planet.Tb_K,3),'%.3f') 'K' ];
        %try
            %load(fullfile([datpath savefile '_pp' num2str]));
            [Planet.Ocean.comp, Planet.FeCore, Planet.Ocean.w_ocean_pct, Planet.Tb_K, ...
            Planet.zb_outerIce_m, Planet.zClath_m, Pb_MPa, PbI_MPa, deltaP, ...
            CMR2mean, Planet.Qmantle_Wm2, Planet.phi_surface, ...
            R_sil_mean_m, R_sil_range_m, rho_sil_mean_kgm3, R_Fe_mean_m, R_Fe_range_m, rho_core_mean_kgm3, ...
            Params.nsteps_clath, Params.nsteps_iceI, nIceIIILithosphere, nIceVLithosphere, nsteps_hydro, ...
            Params.nsteps_mantle, Params.nsteps_core, ...
            R_sil_trade_m, R_Fe_trade_m, rho_sil_trade_kgm3, permeability(1,:,:), ...
            P_Planet_MPa, T_Planet_K, r_Planet_m, phasePlanet, ...
            rho_pPlanet_kgm3, Cp_Planet_JkgK, alpha_Planet_pK, g_Planet_ms2, ...
            phi_Planet_frac, k_S_mPlanet, VP_Planet_kms, VS_Planet_kms, QS_overfgamma_Planet, ...
            Ks_Planet_GPa, Gs_Planet_GPa] ...
                = reloadFromDisk([datpath savestr]);

        z_m = Planet.R_m - r_Planet_m;
        vfluid_kms = VP_Planet_kms;
        nsteps_tot = nsteps_hydro + Params.nsteps_mantle + Params.nsteps_core;
        indSil = nsteps_hydro + 1;
        %catch
            %error(['ERROR: cfg.CALC_NEW=0 but ' savefile ' was not found. Re-run with cfg.CALC_NEW set to 1 to generate the Profile.']);
        %end

        %'nsteps','n_clath','n_iceI','n_ocean'
        n_clath = Params.nsteps_clath;
        n_iceI = Params.nsteps_iceI;
        n_ocean = Params.nsteps_ocean;
        Tb_K = Planet.Tb_K; % This is only for .mat saving

        if ~exist('max_clath_depth')
            max_clath_depth = 1e15;
        end
        P_MPa = P_Planet_MPa;
        T_K = T_Planet_K;
        phase = phasePlanet;
        rho_kgm3 = rho_pPlanet_kgm3;
        Cp = Cp_Planet_JkgK;
        save(fullfile([datpath savefile '_pp']),'P_MPa','Pb_MPa','PbI_MPa','nIceIIILithosphere','T_K','Tb_K','phase','deltaP','wo','rho_kgm3','Cp','nsteps','n_clath','n_iceI','n_ocean','max_clath_depth');
        clear Tb_K
    end
    
%% BEGIN PLOTTING STUFF

lstr_3 = cell(1);
if ~cfg.SKIP_PROFILES % SKIP_PROFILES to be deprecated in favor of a more robust CALC_NEW

    lstr_3 = [math 'T_{b}' nm ': ' num2str(Planet.Tb_K,'%0.1f') ' K'];
    if Planet.FeCore
        % Core/mantle size tradeoff
        set(0, 'CurrentFigure', figs.core);
        clf;hold all
        set(gcf,'Position', [335 133 854 547], 'Name', lbl.corsz)

        plot(R_Fe_trade_m'*1e-3,R_sil_trade_m'*1e-3);
        legend(lstr_3,'Fontsize',lbl.smtext)
        box on
        xlabel([math 'R_{' nm 'Fe}' nm ' (km)'],'Fontsize',lbl.mltext);
        ylabel([math 'R_{' nm 'sil}' nm ' (km)'],'Fontsize',lbl.mltext);
        title(['Fe core ; ' math 'C/MR^2' bnm ' = ' num2str(Planet.Cmeasured) '\pm' num2str(Planet.Cuncertainty) '; ' math ' W ' bnm ' = ' num2str(wo) ' wt%; ' math '\rho_{' bnm 'sil}' bnm ': ' num2str(rho_sil_mean_kgm3,'%0.0f') '; ' math '\rho_{' bnm 'Fe}' bnm ': ' num2str(rho_core_mean_kgm3,'%0.0f')],'Fontsize',lbl.lgtext)

        print(figs.core,Params.savefigformat,fullfile([figpath savebase vcore cfg.xtn]));
    else
        % Mantle radius/density tradeoff with no core
        set(0, 'CurrentFigure', figs.mant);
        clf;hold all
        set(gcf,'Position', [335 133 696 547], 'Name', lbl.mantl)
        plot(rho_sil_trade_kgm3',R_sil_trade_m'*1e-3);
        lstr_3 = [math 'T_{b}' nm ': ' num2str(Planet.Tb_K,'%0.1f') ' K'];
        legend(lstr_3,'Fontsize',lbl.smtext)
        box on
        xlabel([math '\rho_{' nm 'sil}' nm ' (kg m^{-3})'],'Fontsize',lbl.mltext);
        ylabel([math 'R_{' nm 'sil}' nm ' (km)'],'Fontsize',lbl.mltext)
        title(['No Fe core ; ' math 'C/MR^2' bnm ' = ' num2str(Planet.Cmeasured) '\pm' num2str(Planet.Cuncertainty) ';' math ' W ' bnm ' = ' num2str(wo) ' wt%'],'Fontsize',lbl.lgtext)

        print(figs.mant,Params.savefigformat,fullfile([figpath savebase vmant cfg.xtn]));
    end

    if POROUS
        if Params.HOLD
            set(0, 'CurrentFigure', figs.porP);
            set(gcf, 'Name', lbl.porvP)
        else
            set(0, 'CurrentFigure', figs.porP);
            set(gcf, 'Name', [lbl.porvP ' Tb = ' num2str(Planet.Tb_K)])
            clf;
        end
        hold on
%         plot(Pmantle_MPa,[por_in.vp]-[por_out.vp])
%         plot(Pmantle_MPa,por_in.vs-por_out.vs,'--')
        plot(Pmantle_MPa,phi_mantle_frac*100,'LineWidth',cfg.LW_std)
%         plot(Pmantle_MPa,(por_in.den-por_out.den)./por_in.den*100)
%         plot(Pmantle_MPa,(por_in.vp-por_out.vp)./por_in.vp*100)
%         plot(Pmantle_MPa,(por_in.vs-por_out.vs)./por_in.vs*100,'--')
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
            set(0, 'CurrentFigure', figs.porR);
            set(gcf, 'Name', [lbl.porvR ' Tb = ' num2str(Planet.Tb_K)])
            clf;
        end
        hold on
%         plot(Pmantle_MPa,[por_in.vp]-[por_out.vp])
%         plot(Pmantle_MPa,por_in.vs-por_out.vs,'--')
        hl = plot(phi_mantle_frac*100,r_mantle_m*1e-3,'LineWidth',cfg.LW_std);

%         plot((por_in.den-por_out.den)./por_in.den*100,r_mantle_m*1e-3)
%         plot((por_in.vp-por_out.vp)./por_in.vp*100,r_mantle_m*1e-3)
%         plot((por_in.vs-por_out.vs)./por_in.vs*100,r_mantle_m*1e-3,'--')
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
            set(0, 'CurrentFigure', figs.perm);
            set(gcf, 'Name', [lbl.perme ' Tb = ' num2str(Planet.Tb_K)])
            clf;
        end
        hold on
%         plot(Pmantle_MPa,[por_in.vp]-[por_out.vp])
%         plot(Pmantle_MPa,por_in.vs-por_out.vs,'--')
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
    end % POROUS
end % ~cfg.SKIP_PROFILES

    
    if ~cfg.SKIP_PROFILES
        if Seismic.DO_SEISMIC
    %% plot the seismic data and attenuation
        if Params.HOLD
            set(0, 'CurrentFigure', figs.seis);
            set(gcf,'Position', [476   662   857   505],'Name',lbl.seism)
        else
            set(0, 'CurrentFigure', figs.seis);
            set(gcf,'Position', [476   662   857   505],'Name',[lbl.seism ' Tb = ' num2str(Planet.Tb_K)])
            clf;
        end

        hp = subplot(1,3,1);
        if Params.HOLD
            hold on
        end
        plot(VS_Planet_kms',r_Planet_m'*1e-3,...
            VP_Planet_kms',r_Planet_m'*1e-3,'--','LineWidth',cfg.LW_seism)
            set(gcf,'color','w')
            xlabel('Sound Speeds (km s^{-1})','Fontsize',lbl.mltext)
            ylabel([math 'r_{' nm Planet.name '}' nm ' (km)'],'Fontsize',lbl.mltext);
            set(gca,'xlim',[0 1.1*max(VP_Planet_kms)],'ylim',[0 ymax])
            grid on


        hp = subplot(1,3,2);
        if Params.HOLD
            hold on
        end
        if max(P_Planet_MPa)<100
            plot(10*P_Planet_MPa',r_Planet_m'*1e-3,...
                T_Planet_K',r_Planet_m'*1e-3,'--',...
                rho_pPlanet_kgm3',r_Planet_m'*1e-3,'k-.',...
            'LineWidth',cfg.LW_seism);
            xlabel([math 'P' nm ' (bar), ' math 'T' nm ' (K), ' math '\rho' nm ' (kg m^{-3})'],'Fontsize',lbl.mltext);
            xmax = max([max(rho_pPlanet_kgm3) max(10*P_Planet_MPa) max(T_Planet_K)]);
        else
            plot(P_Planet_MPa',r_Planet_m'*1e-3,...
                T_Planet_K',r_Planet_m'*1e-3,'--',...
                rho_pPlanet_kgm3',r_Planet_m'*1e-3,'k-.',...
            'LineWidth',cfg.LW_seism);
            xlabel([math 'P' nm ' (MPa), ' math 'T' nm ' (K), ' math '\rho' nm ' (kg m^{-3})'],'Fontsize',lbl.mltext);
            xmax = max([max(rho_pPlanet_kgm3) max(P_Planet_MPa) max(T_Planet_K)]);
        end
        set(gca,'ylim',[0 ymax],'xlim',[0 xmax]);
        grid on

        subplot(1,3,3)
        if Params.HOLD
            hold on
        end
        hp = plot(QS_overfgamma_Planet',r_Planet_m'*1e-3,'LineWidth',cfg.LW_seism);
        set(gca,'xscale','log','ylim',[0 ymax],'xlim',[10 1e7],'XTick',[10 100 1000 1e4 1e5 1e6 1e7])
        grid on
        xlabel([math 'Q_S/\omega^{\gamma}'],'Fontsize',lbl.mltext)
        set(gcf,'color','w')

        if Params.HOLD
            set(0, 'CurrentFigure', figs.gsks);
            set(gcf,'Position', [476   662   857   505],'Name',lbl.gs_ks)
        else
            set(0, 'CurrentFigure', figs.gsks);
            set(gcf,'Position', [476   662   857   505],'Name',[lbl.gs_ks ' Tb = ' num2str(Planet.Tb_K)])
            clf;
        end
        hold on
        hp = plot(Ks_Planet_GPa',r_Planet_m'*1e-3,...
            Gs_Planet_GPa',r_Planet_m'*1e-3,'--','LineWidth',cfg.LW_seism);
        set(hp,'LineWidth',cfg.LW_seism);
        ylabel([math 'r_{' nm Planet.name '}' nm ' (km)'],'Fontsize',lbl.mltext);
        xlabel([math 'G_S' nm ' and ' math 'K_S' nm ' (GPa)'],'Fontsize',lbl.mltext);
        axis tight

        if ~cfg.HOLD
            
            print(figs.seis,Params.savefigformat,fullfile([figpath thissavestr '_' vseis cfg.xtn]));
            print(figs.gsks,Params.savefigformat,fullfile([figpath thissavestr '_' vgsks cfg.xtn]));
            if POROUS % Saved from earlier -- now we have thissavestr.
                print(figs.porP,cfg.fig_fmt,fullfile([figpath thissavestr '_' vsP   cfg.xtn]));
                print(figs.porR,cfg.fig_fmt,fullfile([figpath thissavestr '_' vsR   cfg.xtn]));
                print(figs.perm,cfg.fig_fmt,fullfile([figpath thissavestr '_' vperm cfg.xtn]));
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
    if Seismic.DO_SEISMIC
        
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

    end
    %% create the legend that describes tb, z_b, and heat flux
    lstr_2 = {};
        lstr_2 = {lstr_2{:} [math 'T_{b}' nm ': ' num2str(Planet.Tb_K,'%0.2f') ' K, ' math 'q_{b}/q_{c}' nm ':'...
            num2str(1e3*Qb,'%0.0f') '/' num2str(1e3*Q_Wm2,'%0.0f') ' mW m^{-2}, ' math 'z_{b}' nm ':' num2str(1e-3*Planet.zb_outerIce_m,'%0.0f') ' km']};

    set(0, 'CurrentFigure', figs.wedg);
    clf;
    subplot(1,1,1); hold on
    input.phase = phasePlanet;
    input.r_m = r_Planet_m;
    PlanetWedge(input);
    box off
    ylabel(Planet.name,'FontSize',24)
    title([math 'T_b' nm ' = ' num2str(Planet.Tb_K) ' K'],'FontSize',lbl.lgtext)
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
    
    if Seismic.DO_SEISMIC
        subplot(2,6,4:6);
    else
        subplot(2,2,2);
    end
    hold on
    line(T_K(1:nsteps_hydro),z_m(1:nsteps_hydro)*1e-3,'Color',Params.colororder(:,1),...
        'LineWidth',cfg.LW_seism,'LineStyle',Params.LineStyle);
    hm = line(T_K(nsteps_hydro),z_m(nsteps_hydro)*1e-3,'Color',Params.colororder(:,1),'Marker','o');
    if strcmp(Params.LineStyle,'-')
        hm.MarkerFaceColor=Params.colororder(:,1);
    end
    set(gca,'ydir','reverse',...
        'xlim',[245 maxScale*max(T_K(1:nsteps_hydro))],'ylim',[0 maxScale*max(Dsil_km)],...
        'FontSize',lbl.mdtext,'YAxisLocation','right');%'XAxisLocation','top');
    set(gcf, 'Name', lbl.panl4);
    
    xlabel('Temperature (K)','Fontsize',lbl.mltext);
    ylabel('Depth (km)','Fontsize',lbl.mltext);
    % for iT=1:nTbs
    %     hline(Dsil_km,Params.colororder(:,1));
    % end
    box on;

    %== Plot sounds speeds in the ice and liquid layers
    if Seismic.DO_SEISMIC
        subplot(2,6,10);
        hold on;
        % clear hp; % putting this here because an error was raised when running
        % Triton models. May 4 2020. I think the problem is the model I'm running
        % doesn't have any liquid. yep.
        hp =  line(vfluid_kms(Params.nsteps_iceI+1:indSil),z_m(Params.nsteps_iceI+1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        maxV = max(vfluid_kms(Params.nsteps_iceI+1:indSil));
        minV = min(vfluid_kms(Params.nsteps_iceI+1:indSil));
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

           subplot(2,6,11);

        hold on;
        line(velsIce.VIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VVt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VVIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
           line(velsIce.Vclatht_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
       ax = gca;
        if Params.HOLD
        
       maxV = max([velsIce.VIt_kms(1:indSil) ...
            velsIce.VIIt_kms(1:indSil) ...
            velsIce.Vclatht_kms(1:indSil)...
            velsIce.VIIIt_kms(1:indSil) ...
            velsIce.VVt_kms(1:indSil) ...
            velsIce.VVIt_kms(1:indSil)]);
       minV = min([velsIce.VIt_kms(1:indSil) ...
         velsIce.Vclatht_kms(1:indSil) ...
            velsIce.VIIt_kms(1:indSil) ...
            velsIce.VIIIt_kms(1:indSil) ...
            velsIce.VVt_kms(1:indSil) ...
            velsIce.VVIt_kms(1:indSil)]);
        maxV = max(maxV);minV = min(minV);

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


        subplot(2,6,12);
        hold on;
        line(velsIce.VIl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VIIIl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        line(velsIce.VVl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
         line(velsIce.Vclathl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_sound);
        
         ax = gca;
        if Params.HOLD
        
        maxV = max([velsIce.VIl_kms(1:indSil) ...
                velsIce.Vclathl_kms(1:indSil) ...
                velsIce.VIIl_kms(1:indSil) ...
                velsIce.VIIIl_kms(1:indSil) ...
                velsIce.VVl_kms(1:indSil) ...
                velsIce.VVIl_kms(1:indSil) ...
                velsIce.Vclathl_kms(1:indSil)]);
        minV = min([velsIce.VIl_kms(1:indSil)  ...
                velsIce.Vclathl_kms(1:indSil) ...
                velsIce.VIIl_kms(1:indSil) ...
                velsIce.VIIIl_kms(1:indSil) ...
                velsIce.VVl_kms(1:indSil) ...
                velsIce.VVIl_kms(1:indSil)]);
        maxV = max(maxV);minV = min(minV);

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
        
    end

    % not sure this is used anymore. flagged for deletion
    % for iT=1:nTbs
    %     if ~isnan(velsIce.VVIl_kms(indSil))
    %         velT = velsIce.VVIt_kms(indSil);
    %         velL = velsIce.VVIl_kms(indSil);
    %     elseif  ~isnan(velsIce.VVl_kms(indSil))
    %         velT = velsIce.VVt_kms(indSil);
    %         velL = velsIce.VVl_kms(indSil);
    %     elseif  ~isnan(velsIce.VIIIl_kms(indSil))
    %         velT = velsIce.VIIIt_kms(indSil);
    %         velL = velsIce.VIIIl_kms(indSil);
    %     elseif  ~isnan(velsIce.VIIl_kms(indSil))
    %         velT = velsIce.VIIt_kms(indSil);
    %         velL = velsIce.VIIl_kms(indSil);
    %     elseif  ~isnan(vfluid_kms(indSil))
    %         velT = vfluid_kms(indSil);
    %         velL = vfluid_kms(indSil);
    %     else
    %         velT = velsIce.VIt_kms(indSil);
    %         velL = velsIce.VIt_kms(indSil);
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
        sigma_Sm_ocean = k_S_mPlanet(Params.nsteps_iceI+1:nsteps_hydro);
        z_km_ocean = z_m(Params.nsteps_iceI+1:nsteps_hydro)*1e-3;
       if wo>0
           if Seismic.DO_SEISMIC
               subplot(3,6,16:18);hold on
           else
               subplot(2,2,4);
           end
        switch Planet.Ocean.comp
            case 'MgSO4'
                line(sigma_Sm_ocean,z_km_ocean,...
                    'Color',Params.colororder(:,1),'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_std);
                mink_val = min(sigma_Sm_ocean);
                maxk_val = max(sigma_Sm_ocean);
            case 'Seawater'
                line(sigma_Sm_ocean,z_km_ocean,'Color',Params.colororder(:,1),...
                    'LineStyle',Params.LineStyle,'LineWidth',cfg.LW_std);
                mink_val = min(sigma_Sm_ocean);
                maxk_val = max(sigma_Sm_ocean);
            case 'NH3'
%                     line(k_S_m(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle',LineStyle,'LineWidth',cfg.LW_std);
                mink_val = NaN;
                maxk_val = NaN;
            case 'NaCl'
                mink_val = NaN;
                maxk_val = NaN;
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
    if Seismic.DO_SEISMIC
        subplot(2,6,[1:3 7:9]);
    else
        subplot(2,2,[1 3]);
    end
    hold on;
    % hold all
    %R2ind = C2mean+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
    R2ind = nsteps_hydro+1+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
    Psil_MPa = diag(P_MPa(R2ind));
    ht= plot(P_MPa(1:nsteps_hydro),rho_kgm3(1:nsteps_hydro),'Color',Params.colororder(:,1),...
        'LineWidth',cfg.LW_std,'LineStyle',Params.LineStyle);
    hm = plot(Psil_MPa,interp1(P_MPa(1:nsteps_hydro),rho_kgm3(1:nsteps_hydro),Psil_MPa),'Color',Params.colororder(:,1));
    if ~Params.HOLD
        hm.MarkerFaceColor=Params.colororder(:,1);
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
    ax.XLim(2) = maxScale*max(Pref_MPa);
    if maxScale*max((rho_kgm3(1:nsteps_hydro)))>ax.YLim(2)
        ax.YLim(2) = maxScale*max(max((rho_kgm3(1:nsteps_hydro))));
    end
    if maxScale*max(Psil_MPa)>ax.XLim(2)
        ax.XLim(2) = maxScale*max(Psil_MPa);
    end
    xlabel('Pressure (MPa)','Fontsize',lbl.mltext)
    ylabel('Density (kg m^{-3})','Fontsize',lbl.mltext)

    print(figs.cond,Params.savefigformat,fullfile([figpath savebase vcond cfg.xtn]));

    %% plot profile with 2 subplots
    else
    if Seismic.DO_SEISMIC
        set(0, 'CurrentFigure', figs.cond);
        maxScale = 1.01;
        if ~Params.HOLD
            clf;
        end
            set(gcf,'Position', [411    52   898    751], 'Name', lbl.panl4)
            subplot(1,2,2);
            hold on

        Dsil_km=(Planet.R_m-R_sil_mean_m(:))*1e-3;

        hp =  line(vfluid_kms(Params.nsteps_iceI+1:indSil),z_m(Params.nsteps_iceI+1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIIl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VIIIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVIl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.VVIt_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',cfg.LW_sound);
        line(velsIce.Vclathl_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',SoundLineWidth);
        line(velsIce.Vclatht_kms(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','-.','LineWidth',SoundLineWidth);

        if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY % MJS 2021-10-31: This appears to not be used.
            line(k_S_mMgSO4p01Planet(1:indSil)*25,z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineStyle','--','LineWidth',cfg.LW_seism);
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
        %     if ~isnan(velsIce.VVIl_kms(indSil))
        %         velT = velsIce.VVIt_kms(indSil);
        %         velL = velsIce.VVIl_kms(indSil);
        %     elseif  ~isnan(velsIce.VVl_kms(indSil))
        %         velT = velsIce.VVt_kms(indSil);
        %         velL = velsIce.VVl_kms(indSil);
        %     elseif  ~isnan(velsIce.VIIIl_kms(indSil))
        %         velT = velsIce.VIIIt_kms(indSil);
        %         velL = velsIce.VIIIl_kms(indSil);
        %     elseif  ~isnan(velsIce.VIIl_kms(indSil))
        %         velT = velsIce.VIIt_kms(indSil);
        %         velL = velsIce.VIIl_kms(indSil);
        %     elseif  ~isnan(vfluid_kms(indSil))
        %         velT = vfluid_kms(indSil);
        %         velL = vfluid_kms(indSil);
        %     else
        %         velT = velsIce.VIt_kms(indSil);
        %         velL = velsIce.VIt_kms(indSil);
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

        line(ax2,T_K(1:indSil),z_m(1:indSil)*1e-3,'Color',Params.colororder(:,1),'LineWidth',cfg.LW_seism);
        hm = line(T_K(indSil),z_m(indSil)*1e-3,'Color',Params.colororder(:,1),'Marker','o');
        hm.MarkerFaceColor=Params.colororder(:,1);

        set(gca,'ydir','reverse','xlim',[245 320],'ylim',[0 maxScale*max(Dsil_km)]);
        xlabel('Temperature (K)','Fontsize',lbl.mltext);
        ylabel('Depth (km)','Fontsize',lbl.mltext);
        % for iT=1:nTbs
        %     hline(Dsil_km,Params.colororder(:,1));
        % end
        box on;

        subplot(1,2,1);
        hold all
        npre = 0; % npre is essentially set to Params.nsteps_mantle - Params.nsteps_mantle much earlier.
        R2ind = nsteps_hydro+1+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
        %R2ind = C2mean+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
        Psil_MPa = diag(P_MPa(R2ind));
        ht=  plot(P_MPa(1:indSil),rho_kgm3(1:indSil),'Color',Params.colororder(:,1),'LineWidth',cfg.LW_seism);
        hm = plot(Psil_MPa,interp1(P_MPa(:),rho_kgm3(:),Psil_MPa),'Color',Params.colororder(:,1));
        hm.MarkerFaceColor=Params.colororder(:,1);

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
        %     line(P_MPa(:),vfluid_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VIl_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VIt_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VIIl_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VIIt_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VIIIl_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VIIIt_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VVl_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VVt_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VVIl_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        %     line(P_MPa,velsIce.VVIt_kms,'Color',Params.colororder(:,1),'Parent',ax(2),'LineStyle','-.','LineWidth',cfg.LW_sound);
        % end
        % set(ax(2),'ylim',[1.5 5])
        % set(ax(2),'xlim',[0 1800])
        % ylabel('Sound Speeds (km s^{-1})')

            print(figs.cond,Params.savefigformat,fullfile([figpath savefile cfg.xtn]));
    end
    end
    %%

    set(0, 'CurrentFigure', figs.grav);
    set(gcf, 'Name', lbl.gravg)
    if ~Params.HOLD
        clf;
    end
    gp = g_Planet_ms2(nsteps_hydro+1);

    subplot(1,2,1);
    if Params.HOLD
        hold on;
    end
    % plot(g_ms2',r_m'*1e-3,gp,R_sil_mean_m*1e-3,'o')
    hl = plot(g_Planet_ms2',r_Planet_m'*1e-3);
    hold on
    hp = plot(gp,R_sil_mean_m*1e-3,'o');
    set(hl,'LineWidth',cfg.LW_seism,'LineStyle',Params.LineStyle,'Color',Params.colororder(:,1));
    set(hp,'Color',Params.colororder(:,1));
    xlabel('g (m s^{-2})','Fontsize',lbl.mltext)
    ylabel('R (km)','Fontsize',lbl.mltext)
    axis tight

    subplot(1,2,2);
    if Params.HOLD
        hold on;
    end
    hl = plot(P_Planet_MPa',r_Planet_m'*1e-3);
    set(hl,'LineWidth',cfg.LW_seism,'LineStyle',Params.LineStyle,'Color',Params.colororder(:,1));

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
    %     ht=  plot(P_MPa,T_K,Params.colororder(:,1));
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
    if Params.CALC_NEW
        endStr = '_new';
    else
        endStr = '_old';
    end
    if cfg.TESTING
        close all
        save(['run_end_' Planet.Ocean.comp endStr])
    end
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
                %disp('WARNING: extrapolating fluid rho, Cp, and alpha, but this is only valid for Seawater above 120 MPa!')
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
            %out = SeaFreeze({P_MPa T_K mo},'NH3'); % MJS 2021-12-03: This doesn't work; SeaFreeze appears to not support anything but pure water as of version 0.9.2.
            disp('WARNING: NH3 composition is not implemented in SeaFreeze, modeling pure water only.')
            out = SeaFreeze({P_MPa T_K},'water1');
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
    if Xin
        if mod(Xin,1)
            d_str = [' &' num2str(Xin,'%0.2f') ];
        else
            d_str = [' &' num2str(Xin,'%0.0f') ];
        end
    else
        d_str = '& -';
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

function writeToDisk(fileName, comp, FeCore, wtPpt, Tb_K, zb_km, zClath_m, Pb_MPa, PbI_MPa, deltaP, ...
    CMR2mean, QfromMantle_Wm2, phiRockMax, RsilMean_m, RsilRange_m, rhoSilMean_kgm3, RFeMean_m, RFeRange_m, rhoCoreMean_kgm3, ...
    nStepsClath, nStepsIceI, nIceIIILitho, nIceVLitho, nStepsHydro, nStepsSil, nStepsCore, ...
    RsilTrade_m, RFeTrade_m, rhoSilTrade_kgm3, permeability, ...
    P_Planet_MPa, T_Planet_K, r_Planet_m, phasePlanet, ...
    rho_pPlanet_kgm3, Cp_Planet_JkgK, alpha_Planet_pK, g_Planet_ms2, ...
    phi_Planet_frac, sig_Planet_Sm, VP_Planet_kms, VS_Planet_kms, QS_Planet, ...
    Ks_Planet_GPa, Gs_Planet_GPa)

    saveFile = fullfile([fileName '.txt']);
    mantCoreFile = fullfile([fileName '_mantleCore.txt']);
    permFile = fullfile([fileName '_mantlePerm.txt']);

    dlmwrite(saveFile, '  nHeadLines = 27', 'delimiter', '');
    dlmwrite(saveFile, ['  Ocean salt = ' comp], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  Iron core = ' num2str(FeCore)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  Salinity(ppt) = ' num2str(wtPpt)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  Tb_K = ' num2str(Tb_K)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  zb_km = ' num2str(zb_km)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  zClath_m = ' num2str(zClath_m)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  Pb_MPa = ' num2str(Pb_MPa)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  PbI_MPa = ' num2str(PbI_MPa)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  deltaP = ' num2str(deltaP)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  CMR2mean = ' num2str(CMR2mean)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  QfromMantle_Wm2 = ' num2str(QfromMantle_Wm2)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  phiRockMax = ' num2str(phiRockMax)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  RsilMean_m = ' num2str(RsilMean_m)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  RsilRange_m = ' num2str(RsilRange_m)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  rhoSilMean_kgm3 = ' num2str(rhoSilMean_kgm3)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  RFeMean_m = ' num2str(RFeMean_m)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  RFeRange_m = ' num2str(RFeRange_m)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  rhoCoreMean_kgm3 = ' num2str(rhoCoreMean_kgm3)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  nStepsClath = ' num2str(nStepsClath)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  nStepsIceI = ' num2str(nStepsIceI)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  nIceIIILitho = ' num2str(nIceIIILitho)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  nIceVLitho = ' num2str(nIceVLitho)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  nStepsHydro = ' num2str(nStepsHydro)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  nStepsSil = ' num2str(nStepsSil)], 'delimiter', '', '-append');
    dlmwrite(saveFile, ['  nStepsCore = ' num2str(nStepsCore)], 'delimiter', '', '-append');
    header = sprintf(['%s' repmat('\t\t%s', 1, 14)],...
    'P (MPa)', 'T (K)', 'r (m)', 'phase ID', 'rho (kg/m3)', ...
    'Cp (J/kg/K)', 'alpha (1/K)', 'g (m/s2)', 'phi (void/solid frac)', ...
    'sigma (S/m)', 'VP (km/s)', 'VS (km/s)', 'QS', 'KS (GPa)', 'GS (GPa)');
    dlmwrite(saveFile, header, 'delimiter', '', '-append');

    saveData = [P_Planet_MPa' T_Planet_K' r_Planet_m' phasePlanet' rho_pPlanet_kgm3' ...
        Cp_Planet_JkgK' alpha_Planet_pK' g_Planet_ms2' phi_Planet_frac' ...
        sig_Planet_Sm' VP_Planet_kms' VS_Planet_kms' QS_Planet' Ks_Planet_GPa' Gs_Planet_GPa'];
    dlmwrite(saveFile,saveData,...
        'delimiter','\t',...
        'precision',18,...
        '-append');
    
    mantCoreHeader = sprintf('%s\t%s\t%s\t', 'RsilTrade (m)', 'RFeTrade (m)', 'rhoSilTrade (kg/m3)');
    dlmwrite(mantCoreFile, mantCoreHeader, 'delimiter', '');
    mantCoreData = [RsilTrade_m' RFeTrade_m' rhoSilTrade_kgm3'];
    dlmwrite(mantCoreFile,mantCoreData,...
        'delimiter','\t',...
        'precision',18,...
        '-append');
    
    permHeader = sprintf('%s\t%s\t%s\t%s\t%s\t', ...
        'permeability1', 'permeability2', 'permeability3', 'permeability4', 'permeability5');
    dlmwrite(permFile, permHeader, 'delimiter', '');
    permData = [permeability(1,:,1)' permeability(1,:,2)' permeability(1,:,3)' ...
         permeability(1,:,4)' permeability(1,:,5)'];
    dlmwrite(permFile,permData,...
        'delimiter','\t',...
        'precision',18,...
        '-append');
end

function [comp, FeCore, wtPct, Tb_K, zb_m, zClath_m, Pb_MPa, PbI_MPa, deltaP, CMR2mean, ...
    QfromMantle_Wm2, phiRockMax, RsilMean_m, RsilRange_m, rhoSilMean_kgm3, RFeMean_m, RFeRange_m, rhoCoreMean_kgm3, ...
    nStepsClath, nStepsIceI, nIceIIILitho, nIceVLitho, nStepsHydro, nStepsSil, nStepsCore, ...
    RsilTrade_m, RFeTrade_m, rhoSilTrade_kgm3, permeability, ...
    P_Planet_MPa, T_Planet_K, r_Planet_m, phasePlanet, ...
    rho_pPlanet_kgm3, Cp_Planet_JkgK, alpha_Planet_pK, g_Planet_ms2, ...
    phi_Planet_frac, sig_Planet_Sm, VP_Planet_kms, VS_Planet_kms, QS_Planet, ...
    Ks_Planet_GPa, Gs_Planet_GPa] ...
        = reloadFromDisk(fileName)

    saveFile = fullfile([fileName '.txt']);
    mantCoreFile = fullfile([fileName '_mantleCore.txt']);
    permFile = fullfile([fileName '_mantlePerm.txt']);
    
    % Parse header from text file
    fReload = fopen(saveFile);
        nHeadLinesStr = split(fgetl(fReload),'=');
        compSplit = split(fgetl(fReload),'=');
        FeCoreStr = split(fgetl(fReload),'=');
        wtPptStr = split(fgetl(fReload),'=');
        Tb_KStr = split(fgetl(fReload),'=');
        zb_kmStr = split(fgetl(fReload),'=');
        zClath_mStr = split(fgetl(fReload),'=');
        Pb_MPaStr = split(fgetl(fReload),'=');
        PbI_MPaStr = split(fgetl(fReload),'=');
        deltaPStr = split(fgetl(fReload),'=');
        CMR2meanStr = split(fgetl(fReload),'=');
        QfromMantle_Wm2Str = split(fgetl(fReload),'=');
        phiRockMaxStr = split(fgetl(fReload),'=');
        RsilMean_mStr = split(fgetl(fReload),'=');
        RsilRange_mStr = split(fgetl(fReload),'=');
        rhoSilMean_kgm3Str = split(fgetl(fReload),'=');
        RFeMean_mStr = split(fgetl(fReload),'=');
        RFeRange_mStr = split(fgetl(fReload),'=');
        rhoCoreMean_kgm3Str = split(fgetl(fReload),'=');
        nStepsClathStr = split(fgetl(fReload),'=');
        nStepsIceIStr = split(fgetl(fReload),'=');
        nIceIIILithoStr = split(fgetl(fReload),'=');
        nIceVLithoStr = split(fgetl(fReload),'=');
        nStepsHydroStr = split(fgetl(fReload),'=');
        nStepsSilStr = split(fgetl(fReload),'=');
        nStepsCoreStr = split(fgetl(fReload),'=');

        nHeadLines = str2double(nHeadLinesStr{2});
        comp = strip(compSplit{2});
        wtPpt = str2double(wtPptStr{2});
        Tb_K = str2double(Tb_KStr{2});
        zb_m = str2double(zb_kmStr{2}) * 1e3; % Value stored as km
        zClath_m = str2double(zClath_mStr{2});
        Pb_MPa = str2double(Pb_MPaStr{2});
        PbI_MPa = str2double(PbI_MPaStr{2});
        deltaP = str2double(deltaPStr{2});
        CMR2mean = str2double(CMR2meanStr{2});
        QfromMantle_Wm2 = str2double(QfromMantle_Wm2Str{2});
        phiRockMax = str2double(phiRockMaxStr{2});
        RsilMean_m = str2double(RsilMean_mStr{2});
        RsilRange_m = str2double(RsilRange_mStr{2});
        rhoSilMean_kgm3 = str2double(rhoSilMean_kgm3Str{2});
        RFeMean_m = str2double(RFeMean_mStr{2});
        RFeRange_m = str2double(RFeRange_mStr{2});
        rhoCoreMean_kgm3 = str2double(rhoCoreMean_kgm3Str{2});
        nStepsClath = str2double(nStepsClathStr{2});
        nStepsIceI = str2double(nStepsIceIStr{2});
        nIceIIILitho = str2double(nIceIIILithoStr{2});
        nIceVLitho = str2double(nIceVLithoStr{2});
        nStepsHydro = str2double(nStepsHydroStr{2});
        nStepsSil = str2double(nStepsSilStr{2});
        nStepsCore = str2double(nStepsCoreStr{2});
    fclose(fReload);
    
    if strcmp(comp, 'Seawater')
        wtPct = wtPpt; % Matlab version uses ppt for Seawater in w_ocean_pct variable
    else
        wtPct = wtPpt/10;
    end

    reloadData = dlmread(saveFile, '', nHeadLines);

    P_Planet_MPa = reloadData(:,1);
    T_Planet_K = reloadData(:,2);
    r_Planet_m = reloadData(:,3);
    phasePlanet = reloadData(:,4);
    rho_pPlanet_kgm3 = reloadData(:,5);
    Cp_Planet_JkgK = reloadData(:,6);
    alpha_Planet_pK = reloadData(:,7);
    g_Planet_ms2 = reloadData(:,8);
    phi_Planet_frac = reloadData(:,9);
    sig_Planet_Sm = reloadData(:,10);
    VP_Planet_kms = reloadData(:,11);
    VS_Planet_kms = reloadData(:,12);
    QS_Planet = reloadData(:,13);
    Ks_Planet_GPa = reloadData(:,14);
    Gs_Planet_GPa = reloadData(:,15);
    
    mantCoreData = dlmread(mantCoreFile, '', 1);
    RsilTrade_m = mantCoreData(:,1);
    if strcmp(strip(FeCoreStr{2}),'True') || strcmp(strip(FeCoreStr{2}),'1')
        FeCore = 1;
        RFeTrade_m = mantCoreData(:,2);
        rhoSilTrade_kgm3 = zeros(size(RsilTrade_m));
    else
        FeCore = 0;
        RFeTrade_m = zeros(size(RsilTrade_m));
        rhoSilTrade_kgm3 = mantCoreData(:,3);
    end
    
    permeability = zeros(nStepsSil,5);
    permData = dlmread(permFile, '', 1);
    permeability(:,1) = permData(:,1);
    permeability(:,2) = permData(:,2);
    permeability(:,3) = permData(:,3);
    permeability(:,4) = permData(:,4);
    permeability(:,5) = permData(:,5);
end