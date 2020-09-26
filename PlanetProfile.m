function PlanetProfile(Planet,Seismic,Params)
%% PlanetProfile
%
% Based on 
% S. Vance, M. Bouffard, M. Choukroun, and C. Sotin. 
% Ganymede's internal structure including thermodynamics of magnesium sulfate oceans in contact with ice. 
% Planetary And Space Science, 96:62-70, 2014.    
% (http://dx.doi.org/10.1016/j.pss.2014.03.011)]
%
% Expanded in
% S. D. Vance, M. P. Panning, S. Staehler, F. Cammarano, B. G. Bills, G. Tobie, S..
% Kamata, S. Kedar, C. Sotin, W. T. Pike, et al. 
% Geophysical investigations of habitability in ice-covered ocean worlds. 
% Journal of Geophysical Research: Planets, Nov 2018.


set(0,'defaultfigurecolor',[1 1 1]) % white background for figures

%% globals for functions in fzero
global M_Planet_kg R_Planet_m extrapflag
extrapflag = 0; % for seawater thermodynamics above 120 MPa
    wo = Planet.Ocean.w_ocean_pct;
if strcmp(Planet.Ocean.comp,'Seawater')
    global swEOS
    swEOS.gsw = swEOS_chooser('gsw305');
elseif strcmp(Planet.Ocean.comp,'NaCl')
    global swEOS
     swEOS.NaCl = readh5spline('NaCl.h5');
     load('NaCl_LBF');
     swEOS.NaCl.sp = sp_NaCl_8GPa;
%    load('NaClLBF');
%    swEOS.NaCl.sp = NaCl_LBF; % up to 4GPa and down to 250 K
elseif strcmp(Planet.Ocean.comp,'NH3')
    %comment out the line below if your have REFPROP installed.
%     error(['NH3 is not currently implemented due to complications with Refprop software.'])
elseif strcmp(Planet.Ocean.comp,'MgSO4')
    conduct_scaling_MgSO4 = (1+4*wo); % empirical scaling of electrical conductivity from 1 bar values compiled in Hand and Chyba 2007
end

thiseos = split(Seismic.mantleEOS,'.tab');
thiseos = char(thiseos(1));
bar2GPa = 1e-4;

M_Planet_kg = Planet.M_kg;
R_Planet_m = Planet.R_m;
Gg = 6.67300e-11; % m3 kg-1 s-2

% implementing a feature to track silicate composition in output files
% because it's getting confusing as we investigate k2, Q, etc.... 
if isfield(Seismic,'mantleEOSname')
    minEOS = Seismic.mantleEOSname;
else
    minEOS = [];
end
savefile = [Planet.name 'Profile_' Planet.Ocean.comp ...
    num2str(round(Planet.Ocean.w_ocean_pct)) 'WtPct_' minEOS];

datpath = strcat(Planet.name,'/');
figpath = strcat(Planet.name,'/figures/');

if ~isfield(Planet,'NoH2O') % backward compatibility--haven't finished implementing water-free worldsc
    Planet.NoH2O =0;
end
%%
if ~Planet.NoH2O
    n_iceI=Params.nsteps_iceI;
    n_ocean = Params.nsteps_ocean;

    if Params.CALC_NEW
    nTbs = length(Planet.Tb_K);

    nsteps = n_iceI+n_ocean;
    [T_K,P_MPa,rho_kgm3] = deal(zeros(nTbs,nsteps));

    phase = zeros(nTbs,nsteps);
    phase(:,1:n_iceI)=1;
    Tb_K = Planet.Tb_K;


    %%
    %ice Ih, ice III, ice V
    %--------------------------------------------------------------------------
    for kt = 1:nTbs %draw four thermal profiles corresponding to four different choices of temperature at the bottom of the Ice I shell
        T_K(kt,1) = Planet.Tsurf_K;
        P_MPa(kt,1) = Planet.Psurf_MPa;

        rho_kgm3(kt,1) = getRhoIce(P_MPa(kt,1),T_K(kt,1),1);
        try
            Pb_MPa(kt) = getPfreeze(Tb_K(kt),wo,Planet.Ocean.comp);
            deltaP = Pb_MPa(kt)/(n_iceI-1);            

            for il=2:n_iceI % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                P_MPa(kt,il) = P_MPa(kt,il-1) + deltaP;   
                T_K(kt,il) = (Planet.Tb_K(kt).^(P_MPa(kt,il)./Pb_MPa(kt))).*(Planet.Tsurf_K.^(1-P_MPa(kt,il)./Pb_MPa(kt)));
                rho_kgm3(kt,il) = getRhoIce(P_MPa(kt,il),T_K(kt,il),1); 
            end
            nIceIIILithosphere=0;
            PbI_MPa = Pb_MPa;
        catch
            disp('PlanetProfile failed to get Pb! Here''s why:')
            disp(lasterr)            
            disp('Maybe it''s okay? If execution stopped, then probably not. Try looking where getPfreeze is called.')
            
            if isfield(Params,'BOTTOM_ICEIII') && Params.BOTTOM_ICEIII % this will elicit an error if one has set the temperature too low but one hasn't specified ice III or V under the ice I
                disp('Adding ice III to the bottom of ice Ih. Make sure the ocean salinity is high enough that doing this makes sense')
                nIceIIILithosphere=5;
                phase(kt,n_iceI-nIceIIILithosphere:n_iceI)=3; 

                PbI_MPa(kt) = 210; % the Ih-III transition is essentially fixed, within a few MPa    
                deltaP = PbI_MPa(kt)/(n_iceI-5-1);  % save five entries at the bottom for ice III

                for il=2:n_iceI-nIceIIILithosphere % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                    P_MPa(kt,il) = P_MPa(kt,il-1) + deltaP;   
                    T_K(kt,il) = (Planet.Tb_K(kt).^(P_MPa(kt,il)./PbI_MPa(kt))).*(Planet.Tsurf_K.^(1-P_MPa(kt,il)./PbI_MPa(kt)));
                    rho_kgm3(kt,il) = getRhoIce(P_MPa(kt,il),T_K(kt,il),1); 
                end

                TbIII = Tb_K(kt)+1;
                Pb_MPa(kt) = getPfreezeIII(TbIII,wo,Planet.Ocean.comp);% fix the thickness of the ice III layer based on T>Tb by 2K
                deltaP = (Pb_MPa(kt)-PbI_MPanIceIIILithosphere);  %  five entries at the bottom for ice III
                for il = n_iceI-nIceIIILithosphere+1:n_iceI
                    P_MPa(kt,il) = P_MPa(kt,il-1) + deltaP;   
                    T_K(kt,il) = (TbIII.^(P_MPa(kt,il)./Pb_MPa(kt))).*(TbIII(kt).^(1-P_MPa(kt,il)./Pb_MPa(kt)));
                    rho_kgm3(kt,il) = getRhoIce(P_MPa(kt,il),T_K(kt,il),3); 
                end
            elseif isfield(Params,'BOTTOM_ICEV') && Params.BOTTOM_ICEV
                 disp('Adding ice V and III to the bottom of ice Ih. Make sure the ocean salinity is high enough that doing this makes sense')
                nIceIIILithosphere=5;
                nIceVLithosphere=5;
                phase(kt,n_iceI-nIceVLithosphere-nIceIIILithosphere:n_iceI-nIceVLithosphere)=3; 
                phase(kt,n_iceI-nIceVLithosphere:n_iceI)=5; 

                PbI_MPa(kt) = 210; % the Ih-V transition is essentially fixed, within a few MPa    
                deltaP = PbI_MPa(kt)/(n_iceI-5-1);  % save five entries at the bottom for ice V

                %ice Ih
                for il=2:n_iceI-nIceVLithosphere % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
                    P_MPa(kt,il) = P_MPa(kt,il-1) + deltaP;   
                    T_K(kt,il) = (Planet.Tb_K(kt).^(P_MPa(kt,il)./PbI_MPa(kt))).*(Planet.Tsurf_K.^(1-P_MPa(kt,il)./PbI_MPa(kt)));
                    rho_kgm3(kt,il) = getRhoIce(P_MPa(kt,il),T_K(kt,il),1); 
                end

                %ice III
                TbIII = Tb_K(kt)+1;
                Pb_MPa(kt) = getPfreezeIII(TbIII,wo,Planet.Ocean.comp);% fix the thickness of the ice III layer based on T>Tb by 2K
                deltaP = (Pb_MPa(kt)-PbI_MPa)/(nIceIIILithosphere);  %  five entries at the bottom for ice III
                for il = n_iceI-nIceIIILithosphere+1:n_iceI
                    P_MPa(kt,il) = P_MPa(kt,il-1) + deltaP;   
                    T_K(kt,il) = (TbIII.^(P_MPa(kt,il)./Pb_MPa(kt))).*(TbIII(kt).^(1-P_MPa(kt,il)./Pb_MPa(kt)));
                    rho_kgm3(kt,il) = getRhoIce(P_MPa(kt,il),T_K(kt,il),3); 
                end
                
                %ice V
                TbV = Tb_K(kt)+1;
                Pb_MPa(kt) = getPfreezeV(TbV,wo,Planet.Ocean.comp);% fix the thickness of the ice V layer based on T>Tb by 2K
                deltaP = (Pb_MPa(kt)-PbI_MPa)/(nIceVLithosphere);  %  five entries at the bottom for ice V
                for il = n_iceI-nIceVLithosphere+1:n_iceI
                    P_MPa(kt,il) = P_MPa(kt,il-1) + deltaP;   
                    T_K(kt,il) = (TbV.^(P_MPa(kt,il)./Pb_MPa(kt))).*(TbV(kt).^(1-P_MPa(kt,il)./Pb_MPa(kt)));
                    rho_kgm3(kt,il) = getRhoIce(P_MPa(kt,il),T_K(kt,il),3); 
                end
            end
        end



    %OCEAN + ICE III/V/VI SHELL
    %--------------------------------------------------------------------------
      deltaP = (Params.Pseafloor_MPa-Pb_MPa(kt))/n_ocean; %
      if deltaP<=0 % prevent an embarassing error that can occur when adapting a file for a small object to that of a larger one.
          error('negative increment of pressure while computing ocean profile. Be sure Params.Pseafloor_MPa is greater than the likely pressure at the silicate interface. Currently it''s less than the pressure at the base of the ice I layer.')
      end
      for il =1:n_ocean
        ill = il+n_iceI;
        disp(['kt: ' num2str(kt) '; il: ' num2str(il) '; P_MPa: ' num2str(round(P_MPa(kt,ill-1)))]);
        if il==178
            break_here = 1;
        end
        P_MPa(kt,ill) = Pb_MPa(kt) + il*deltaP;
        
        if il==1 % establish the phase vector
            phase(kt,ill) = getIcePhase(P_MPa(kt,ill),T_K(kt,ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
        elseif phase(kt,ill-1)~=6
            phase(kt,ill) = getIcePhase(P_MPa(kt,ill),T_K(kt,ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
            if phase(kt,ill-1)==3 && phase(kt,ill)==0 % fix to instabilities in the phase lookup for ice III
                disp('Fixing apparent instability in the EOS?ice III is forming above the ocean where it wasn''t specified')
                phase(kt,ill)=3;
            end
            if phase(kt,ill-1)==0 && phase(kt,ill)==1 % fix to instabilities in the phase lookup for ice III
                disp('Fixing apparent instability in the EOS?ice Ih is forming under the ocean')
                phase(kt,ill)=0;
            end
        else
            phase(kt,ill) = 6;
        end


        if phase(kt,ill) == 0 %if ocean
          [rho_ocean,Cp,alpha_o]= fluidEOS(P_MPa(kt,ill),T_K(kt,ill-1),wo,Planet.Ocean.comp);
          if alpha_o<=0
              disp('Ocean alpha at ice interface is less than zero. adjusting temperature upward.')
              disp('This means there''s a conductive layer at the interface with thickness inversely proportional to the heat flow.')
              disp('The thickness is likely less than a few 100 m. See Melosh et al. 2004.')
              Tnew = fzero(@(T) alphaAdjust(P_MPa(kt,ill),T,wo,Planet.Ocean.comp),273);
              disp(['The new temperature is ' num2str(Tnew) ', or delta T = ' num2str(Tnew-T_K(kt,ill-1)) '.'])
              T_K(kt,ill)=Tnew;
          else
              T_K(kt,ill) = T_K(kt,ill-1)+ alpha_o*T_K(kt,ill-1)./(Cp)./(rho_ocean)*deltaP*1e6; %adiabatic gradient in ocean; this introduces an error, in that we are using the temperature from the previous step
          end
          rho_kgm3(kt,ill) = fluidEOS(P_MPa(kt,ill),T_K(kt,ill),wo,Planet.Ocean.comp);    % THIS IS REDUNDANT TO THE CALCULATION OF RHO_OCEAN ABOVE
        else %if ice
            % allow for stable dense fluid under high pressure ice --> this
            % was added for the callisto study. Commented out currently
            rhoice = getRhoIce(P_MPa(kt,ill),T_K(kt,ill-1),phase(kt,ill));
%             rhoocean = fluidEOS(P_MPa(kt,ill),T_K(kt,ill-1),wo,Planet.Ocean.comp);
%             if rhoocean>=rhoice
%                 phase(kt,ill)=0;
%             end
            if ~isfield(Planet.Ocean,'fnTfreeze_K')
                if il==6
                    x = 1;
                end
                T_K(kt,ill) = getTfreeze(P_MPa(kt,ill),wo,Planet.Ocean.comp,T_K(kt,ill-1)); %find the temperature on the liquidus line corresponding to P(k,i+nIceI); should really use a conductive profile here, but that would seem to naturally bring us back to the liquidus. Steve wants to confirm.
            else
                T_K(kt,ill) = Planet.Ocean.fnTfreeze_K(P_MPa(kt,ill),wo);
            end
            if T_K(kt,ill)<T_K(kt,ill-1)
                T_K(kt,ill)=T_K(kt,ill-1);
            end
            if phase(kt,ill)==0
                rho_kgm3(kt,ill) = rhoocean;
            else
                rho_kgm3(kt,ill) = getRhoIce(P_MPa(kt,ill),T_K(kt,ill),phase(kt,ill));
            end
        end
    end
    rho_kgm3(kt,1) = rho_kgm3(kt,2); %continuity
        save(strcat(datpath,savefile),'P_MPa','Pb_MPa','PbI_MPa','nIceIIILithosphere','T_K','Tb_K','phase','deltaP','wo','nTbs','rho_kgm3','rho_ocean','Cp','alpha_o','nsteps'); % save the progress at each step
    end  
    else
        load(strcat(datpath,savefile));
    end
    %% Save in a file the densities of adiabats corresponding to different ocean concentrations
    str_ref_densities = ['ref_densities' Planet.name '_' Planet.Ocean.comp '.mat'];
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
        save([datpath str_ref_densities],'rho_ref_kgm3','Pref_MPa','Tref_K')
    else
        load([datpath str_ref_densities]);
    end

    %%%%%%%%%%%%%%%%%%%%%
    % convert to depth
    %%%%%%%%%%%%%%%%%%%
    %% calculate gravity in each layer instead of assuming surface gravity applies.
    % allocate variables
        z_m = zeros(nTbs,nsteps);
        g_ms2 = z_m;
        [M_above_kg,M_below_kg] = deal(rho_kgm3); % mass above and below the silicate interface
        M_above_kg(:,1) = 0;
        M_below_kg(:,1) = M_Planet_kg;
        [z_m,r_m] = deal(zeros(nTbs,nsteps)); 
        r_m(:,1) = Planet.R_m;
        g_ms2(1:nTbs,1) = Gg*M_Planet_kg/R_Planet_m^2;

    for kt = 1:nTbs
        deltaP = Pb_MPa(kt)/n_iceI;
        for il = 2:n_iceI
            % calculate depth
            z_m(kt,il) = z_m(kt,il-1)+ (P_MPa(kt,il)-P_MPa(kt,il-1))*1e6/g_ms2(kt,il-1)/rho_kgm3(kt,il-1);
            % convert to radius
            r_m(kt,il) = Planet.R_m-z_m(kt,il); 

            % determine local gravity
            M_above_kg(kt,il) = M_above_kg(kt,il-1) + 4/3*pi*(r_m(kt,il-1)^3-r_m(kt,il)^3)*rho_kgm3(kt,il);
            M_below_kg(kt,il) = M_Planet_kg-M_above_kg(kt,il);
            g_ms2(kt,il) = Gg*M_below_kg(kt,il)/r_m(kt,il)^2;
        end
        Zb2(kt)=z_m(kt,il);  
        deltaP = (Params.Pseafloor_MPa-Pb_MPa(kt))/n_ocean; %
        for il = 1:n_ocean
            ill = il+n_iceI;
            %calculate depth
            z_m(kt,ill) = z_m(kt,il-1+n_iceI)+ deltaP*1e6/g_ms2(kt,il-1+n_iceI)/rho_kgm3(kt,ill); % using the previous gravity step, since we haven't calculated the present step.  this introduces an error
            % convert to radius
            r_m(kt,ill) = R_Planet_m-z_m(kt,ill); 

            % determine local gravity
            M_above_kg(kt,ill) = M_above_kg(kt,il-1+n_iceI)+4/3*pi*(r_m(kt,il-1+n_iceI)^3-r_m(kt,ill)^3)*rho_kgm3(kt,ill);
            M_below_kg(kt,ill) = M_Planet_kg-M_above_kg(kt,ill);
            g_ms2(kt,ill) = Gg*M_below_kg(kt,ill)/r_m(kt,ill)^2;
        end
    end

    %% compute conductive heat through the ice I layer
    D_conductivityIh = 632; % W m-1; Andersson et al. 2005 (For comparison, Mckinnon 2006 uses a value of 621 from Slack 1980)
    for kt = 1:nTbs
        Qb(kt) = D_conductivityIh*log(Planet.Tb_K(kt)/Planet.Tsurf_K)/Zb2(kt);
    end
end

%% Deeper Interior
    MR2 = M_Planet_kg*R_Planet_m^2;
    dz_m = gradient(z_m);
    nz = length(dz_m);
    nR = nz;
%     nR = round(3/4*nz);
    npre = nz-nR; % confine the search for structure to the lower part
    C_H2O = sum((8/3*pi*rho_kgm3(:,1:npre).*r_m(:,1:npre).^4.*dz_m(:,1:npre))')'; %calculate C_H2O to the beginning of the depth where we start probing for the beginning of the silicate layer

    
    % without a core
if ~Planet.FeCore
    C1 = zeros(nTbs,nR);
    for iz = 1:nR
        C_H2O(:,iz+1) = C_H2O(:,iz)+(8/3*pi*rho_kgm3(:,iz+npre).*r_m(:,iz+npre).^4.*dz_m(:,iz+npre));
        R_sil_m(:,iz) = r_m(:,iz+npre);
        rho_sil_kgm3(:,iz) = 3/4/pi*(M_Planet_kg-M_above_kg(:,iz+npre))./power(R_sil_m(:,iz),3);
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
   
%     dzOcean_m(~dindsVI & ~dindsV) = R_Planet_m - R_sil_mean_m(~dindsVI &
%     ~dindsV)-zI_m(~dindsVI & ~dindsV); not sure this is right. commenting
%     out on March 8 2018

    
   % Plot the results
if ~Params.NOPLOTS
    figure(2295);clf;hold all
    for iT=1:nTbs
        plot(rho_sil_kgm3(iT,C2inds{iT})',R_sil_m(iT,C2inds{iT})'*1e-3);
    end
    lstr_3 = {};
    for iT = 1:nTbs
        lstr_3 = {lstr_3{:} ['T_{b}:' num2str(Planet.Tb_K(iT),'%0.1f') 'K']};
    end
    legend(lstr_3)
    box on
    xlabel('\rho_{sil} (kg m^{-3})');
    ylabel('R_{sil} (km)')
    title(['No Fe core ; C/MR2=' num2str(Planet.Cmeasured) '\pm' num2str(Planet.Cuncertainty) '; W =' num2str(wo) ' Wt%'])
end
% =====
else % WITH A CORE
    rho_Fe_kgm3 = Planet.rhoFe*Planet.rhoFeS/(Planet.xFeS*(Planet.rhoFe-Planet.rhoFeS)+Planet.rhoFeS); 
    C2 = zeros(nTbs,nR);
    M_iron_kg = zeros(nTbs,nR);
    R_Fe_m = M_iron_kg;
    for iz = 1:nR
        C_H2O(:,iz+1) = C_H2O(:,iz)+(8/3*pi*rho_kgm3(:,iz+npre).*r_m(:,iz+npre).^4.*dz_m(:,iz+npre));    
        R_sil_m(:,iz) = r_m(:,iz+npre);
        for iT = 1:nTbs
            [C2(iT,iz),R_Fe_m(iT,iz)] = CoreSize(Planet.rho_sil_withcore_kgm3,rho_Fe_kgm3,C_H2O(iT,iz+1),M_above_kg(iT,iz+npre),R_sil_m(iT,iz));
        end
    end
if ~Params.NOPLOTS    
    figure(2296);clf;hold all
end
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
if ~Params.NOPLOTS
        plot(R_Fe_m(iT,C2inds{iT})'*1e-3,R_sil_m(iT,C2inds{iT})'*1e-3);
end
    end
if ~Params.NOPLOTS
    lstr_3 = {};
    for iT = 1:nTbs
        lstr_3 = {lstr_3{:} ['T_{b}:' num2str(Planet.Tb_K(iT),'%0.1f') 'K']};
    end

    legend(lstr_3)
    box on
    xlabel('R_{Fe} (km)');
    ylabel('R_{sil} (km)')
    title(['Fe core ; C/MR2=' num2str(Planet.Cmeasured) '\pm' num2str(Planet.Cuncertainty) '; W =' num2str(wo) ' Wt%; \rho_{sil}: ' num2str(Planet.rho_sil_withcore_kgm3,'%0.0f') '; \rho_{Fe}: ' num2str(rho_Fe_kgm3,'%0.0f')])
end
end


%% Display the depths to the various layers
%allocate
    zI_m = Zb2;
    zIII_m = zeros(1,nTbs);
    zV_m = zeros(1,nTbs);
    zVI_m = zeros(1,nTbs);
    indsV = zV_m;
    indsVI = zV_m;
    indsIII = zV_m;

    % figure out the indices for the tops of the different ice layers
for iT = 1:nTbs
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
end
[dzOcean_m,dzIII_m,dzV_m,dzVI_m] = deal(zeros(1,nTbs));
%find the radii at the tops of the different layers
RIII_m =R_Planet_m-zIII_m;
RV_m = R_Planet_m-zV_m;
RVI_m = R_Planet_m-zVI_m;

% find the thicknesses of the layers.  keep adjusting the thicknesses
% accordingly
dindsIII = RIII_m>R_sil_mean_m & zIII_m>0;
dzOcean_m(dindsIII) = zIII_m(dindsIII)-zI_m(dindsIII);

dindsV = RV_m>R_sil_mean_m & zV_m>0;
dzOcean_m(~dindsIII & dindsV) = zV_m(~dindsIII & dindsV)-zI_m(~dindsIII & dindsV); 
dzIII_m(dindsIII & dindsV) = zV_m(dindsV & dindsIII)- zIII_m(dindsV & dindsIII);

dindsVI = RVI_m>R_sil_mean_m & zVI_m>0;
dzOcean_m(dindsVI & ~dindsV) = zVI_m(dindsVI & ~dindsV) - zI_m(dindsVI & ~dindsV);

    % Display Results For the Calculation of Mantle and Core Depths
    % use R_sil_mean_m to calculate thickness of ice VI
    for iT = 1:nTbs
        if ~isempty(dindsVI(iT)) && dindsVI(iT)>0
            dzVI_m(iT) = R_Planet_m-R_sil_mean_m(iT) - zVI_m(iT);
        elseif ~isempty(dindsV(iT)) && dindsV(iT)>0 
            dzVI_m(iT) = 0;
            dzV_m(iT) = RV_m(iT)-R_sil_mean_m(iT) - zVI_m(iT)*dindsVI(iT);
        elseif ~isempty(dindsIII(iT)) && dindsIII(iT)>0
            dzV_m(iT) = 0;
            dzVI_m(iT) = 0;
            dzIII_m(iT) = R_Planet_m-R_sil_mean_m(iT) - zVI_m(iT);
%         elseif ~isempty(dindsII(iT)) && dindsII(iT)>0
%             dzV_m(iT) = 0;
%             dzVI_m(iT) = 0;
%             dzIII_m(iT) = R_Planet_m-R_sil_mean_m(iT) - zVI_m(iT);            
%         end
        end
    end

    dzOcean_m(~dindsVI & ~dindsV) = R_Planet_m - R_sil_mean_m(~dindsVI & ~dindsV)-zI_m(~dindsVI & ~dindsV);

    zTotal_m = zI_m+dzOcean_m+dzIII_m+dzV_m+dzVI_m;

    % this should be elaborated upon to compute the actual mass by
    % integrating the masses of the concentric shells and removing the
    % effect of the salt. This can be easily done using eqst. It will be
    % interesting to compare the mass of the pure water ice and ocean layer
    % with that of the salty layer
    % less precise:
%     MH2O = 4/3*pi*(R_Planet_m.^3 - R_sil_mean_m.^3).*mean(rho_kgm3(:,C2mean));
    % more precise:
    WtH2O = M_H2O_mean_kg/M_Planet_kg; % this is the predicted water in the planet
    dWtH2O = Planet.XH2O-WtH2O; 
    
if Planet.FeCore
    Mcore = 4/3*pi*R_Fe_mean_m.^3.*rho_Fe_kgm3;
    Wtcore = Mcore/M_Planet_kg;
    dWtcore = Planet.xFe_core-Wtcore;
    XS = 32.065/87.91*Planet.xFeS; %mass fraction of S in the core
    Msulfur = Mcore.*XS; % mass of sulfur in the core
    WtS = Msulfur/M_Planet_kg; % as a fraction of Europa's mass
    dWtS = Planet.xFeS_meteoritic*32.065/87.91-WtS;
end
    
    disp(['Tb:                    ' num2str(Planet.Tb_K,'\t%0.2f')])
    disp(['z(km) ice I:           ' num2str(zI_m*1e-3,'\t%0.1f')])    
    disp(['z(km) ice III:         ' num2str(dindsIII.*zIII_m*1e-3,'\t%0.0f')])    
    disp(['z(km) ice V:           ' num2str(dindsV.*zV_m*1e-3,'\t%0.0f')])    
    disp(['z(km) ice VI:          ' num2str(dindsVI.*zVI_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) Ocean:          ' num2str(dzOcean_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) ice III:        ' num2str(dzIII_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) ice V:          ' num2str(dzV_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) ice VI:         ' num2str(dzVI_m*1e-3,'\t%0.0f')])  
    disp(['dz(km) ice V + ice VI: ' num2str((dzV_m+dzVI_m)*1e-3,'\t%0.0f')])  
    if isfield(Planet,'XH2O')
    disp('')
    disp(['Wt Pct of H2O:         ' num2str(100*WtH2O,'\t%0.2f')]);
    disp(['Chondritic input H2O is extra     ' num2str(100*dWtH2O,'\t%0.1f') ' % of ' Planet.name '''s mass'])
    disp('')
    if Planet.FeCore
        disp(['Wt Pct of the core:    ' num2str(100*Wtcore,'\t%0.2f')]);
        disp(['Chondritic input core matter is extra ' num2str(100*dWtcore,'\t%0.1f') ' % of ' Planet.name '''s mass'])
        disp('')
        disp(['Wt Pct S as FeS in core:' num2str(100*XS,'\t%0.2f')]);
        disp(['Wt Pct core S in Europa:' num2str(100*WtS,'\t%0.2f')]);
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


%% Iterative Calculation for Adding Details to the Seismic Profile 
% consider, if appropriate, a convective ice shell of the same thickness
% this introduces an error in the gravity profile and thus the moment of inertia, since the overlying mass will be
% less
for iT = 1:nTbs
    % Ice I
    [Q_Wm2(iT),deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),nu(iT),CONVECTION_FLAG_I(iT)]=...
        ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K(iT),PbI_MPa(iT)/2,Zb2(iT),g_ms2(iT,1),2);
    if ~isfield(Params,'NO_ICEI_CONVECTION') || Params.NO_ICEI_CONVECTION == false
        if CONVECTION_FLAG_I(iT) || (isfield(Params,'FORCE_ICEI_CONVECTION') && Params.FORCE_ICEI_CONVECTION == true)
            %conductive upper layer
            nConvectIce=n_iceI-nIceIIILithosphere-1;
            inds_eTBL = find(z_m(iT,1:nConvectIce)<=eTBL_m(iT));
            Pterm = P_MPa(iT,inds_eTBL)./P_MPa(iT,inds_eTBL(end));
            T_K(iT,inds_eTBL) = (Tc(iT).^(Pterm)).*(Planet.Tsurf_K.^(1-Pterm));
            %convective region
            for iconv = inds_eTBL(end)+1:nConvectIce
    %             rho = 1000./getVspChoukroun2010(P_MPa(iT,iconv),T_K(iT,iconv-1),2);
                rho = getRhoIce(P_MPa(iT,iconv),T_K(iT,iconv-1),1);
                try
                [Cp,alpha_K] = getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),1);
                catch
                    x =1;
                end
    %             aK = 1.56e-4;
                T_K(iT,iconv) = T_K(iT,iconv-1)+alpha_K*T_K(iT,iconv)/Cp/rho*deltaP*1e6;
            end
           % conductive lower layer
           inds_deltaTBL = find(z_m(iT,1:nConvectIce)>=z_m(iT,nConvectIce)-deltaTBL_m(iT));
           T_K(iT,inds_deltaTBL) = (Planet.Tb_K(iT).^(P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT))).*(T_K(iT,inds_deltaTBL(1)-1).^(1-P_MPa(iT,inds_deltaTBL)./PbI_MPa(iT)));   
           rho_kgm3(iT,inds_deltaTBL) = getRhoIce(P_MPa(iT,inds_deltaTBL),T_K(iT,inds_deltaTBL),1); 

           if find(phase(iT,:)>1)
    %         indVI = find(phase(iT,:)==6);   
    %         Ttop = T_K(iT,indsVI(iT));
    %         Tbottom = zVI_m
    %         [Q_Wm2(iT),deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),CONVECTION_FLAG_I(iT)]=...
    %         ConvectionHPIceKalousova2018(Ttop,,PbI_MPa(iT)/2,Zb2(iT),g_ms2(iT,1),2);
           end

           if Planet.EQUIL_Q
               Planet.Qmantle_Wm2(iT) = Q_Wm2(iT);
           else
               Planet.Qmantle_Wm2(iT) = Planet.Qmantle_Wm2(1);
           end
        else
           if Planet.EQUIL_Q
               Planet.Qmantle_Wm2(iT) = Qb(iT);
           else
               Planet.Qmantle_Wm2(iT) = Planet.Qmantle_Wm2(1);
           end

        end
    else
        if Planet.EQUIL_Q
           Planet.Qmantle_Wm2(iT) = Qb(iT);
        else
           Planet.Qmantle_Wm2(iT) = Planet.Qmantle_Wm2(1);
        end
        Q_Wm2(iT) = Planet.Qmantle_Wm2(1); 
    end
end

indSil=zeros(1,nTbs);
for iT = 1:nTbs % determine where the silicate should start
    indSil(iT) = find(R_sil_m(iT,:)==R_sil_mean_m(iT));
end
if nTbs>1 % this makes sure that the mantle length matrices are all the same length.   
    nsteps_mantle = [Params.nsteps_mantle Params.nsteps_mantle*ones(1,nTbs-1)+(indSil(1)-indSil(2:end))];% 
    if find(nsteps_mantle<0)
        error('Problems indexing in the "mantle" layer. Consider increasing nsteps_mantle.');
    end
else
    nsteps_mantle = Params.nsteps_mantle;
end

mantle = loadMantleEOS(Seismic.mantleEOS);
% mprops = load(Seismic.mantleEOS);
if isfield(Seismic,'mantleEOS_dry')
    mantleDry = loadMantleEOS(Seismic.mantleEOS_dry);
end   
    

for iT = 1:nTbs
    mean_rho_ocean(iT) = mean(rho_kgm3(iT,n_iceI+1:indSil(iT)));
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
    POROUS = isfield(Planet,'POROUS_ROCK') && Planet.POROUS_ROCK;
for iT = 1:nTbs
    [Pmantle_MPa,rho_mantle_kgm3,M_above_mantle,VP_mantle_kms,VS_mantle_kms,g_ms2_sil] = deal(zeros(1,nsteps_mantle(iT)));
    Tmantle_K = T_K(iT,indSil(iT));
    r_mantle_m = linspace(R_sil_mean_m(iT),R_Fe_mean_m(iT),nsteps_mantle(iT));

    g_ms2_sil = g_ms2(iT,C2mean(iT))*ones(1,nsteps_mantle(iT));
    MantleHeat = Planet.Qmantle_Wm2(iT)*4*pi*R_Planet_m^2+Planet.QHmantle;
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
%             g_ms2_sil(ij)=(Gg*(M_Planet_kg-M_above_mantle(ij-1))/r_mantle_m(ij-1)^2);
%         else
%             g_ms2_sil(ij) = g_ms2_sil(ij-1);
%         end
        Pmantle_MPa(ij) = Pmantle_MPa(ij-1)+1e-6*(rho_mantle_kgm3(ij-1))*g_ms2_sil(ij)*(r_mantle_m(ij-1)-r_mantle_m(ij));
        M_above_mantle(ij) = M_above_mantle(ij-1)+4/3*pi*(r_mantle_m(ij-1)^3-r_mantle_m(ij)^3)*rho_mantle_kgm3(ij-1);
        
        if isfield(Seismic,'mantleEOS_dry') & rho_mantle_kgm3(ij-1)>rhofn_kgm3(Pmantle_MPa(ij),Tmantle_K(ij))
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
%             g_ms2_sil(ij)=(Gg*(M_Planet_kg-M_above_mantle(ij-1))/r_mantle_m(ij-1)^2);
%         else
%             g_ms2_sil(ij) = g_ms2_sil(ij-1);
%         end
%         Pmantle_MPa(ij) = Pmantle_MPa(ij-1)+1e-6*(rho_mantle_kgm3(ij-1))*g_ms2_sil(ij)*(r_mantle_m(ij-1)-r_mantle_m(ij));
          M_above_mantle(ij) = M_above_mantle(ij-1)+4/3*pi*(r_mantle_m(ij-1)^3-r_mantle_m(ij)^3)*rho_mantle_kgm3(ij-1);
        end
    end
    
    mtest = find(M_above_mantle>M_Planet_kg);
    if mtest
        disp(['Exceeded planet mass at mantle radius of ' num2str(1e-3*r_mantle_m(mtest(1)),'%0.0f') ' km']);
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
        r_core_m(iT,:) = linspace(R_Fe_mean_m(iT),0,Params.nsteps_core);
        g_ms2_core(iT,:)=(Gg*(M_Planet_kg-M_above_core(iT))/r_core_m(iT)^2)*ones(1,Params.nsteps_core);

        Tcore_K(iT,:) = linspace(Tmantle_K(end),1.01*Tmantle_K(end),Params.nsteps_core);
        Pcore_MPa(iT,1) = Pmantle_MPa(nsteps_mantle(iT)); 
        if isfield(Seismic,'coreEOS')
            core = loadMantleEOS(Seismic.coreEOS);
            rhocore_kgm3(iT,1) = core.rho_fn(Pcore_MPa(iT,1),Tcore_K(iT,1));
            
            for ij = 2:Params.nsteps_core
        %          g_ms2_core(iT,ij)=(Gg*(M_Planet_kg-M_above_core(iT,ij-1))/r_core_m(iT,ij-1)^2);
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
        %          g_ms2_core(iT,ij)=(Gg*(M_Planet_kg-M_above_core(iT,ij-1))/r_core_m(iT,ij-1)^2);
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
            if ~endsWith(mfnames{im},'_fn') & ~(strcmp(mfnames{im},'p') | strcmp(mfnames{im},'t'))
                interior(iT).(['fluid_' mfnames{im}]) = mfluids.(mfnames{im+1})(Pmantle_MPa,Tmantle_K);
            end
        end
        for im = 1:length(mpnames) % count only the names, not the interpolating functions or p or t
            if ~endsWith(mpnames{im},'_fn') & ~(strcmp(mpnames{im},'p') | strcmp(mpnames{im},'t'))
                interior(iT).(['wt_' mpnames{im}]) = mphases.(mpnames{im+1})(Pmantle_MPa,Tmantle_K); %wt pct of each constituent
            end
        end
        masscheck = 0;
        for im = 1:length(mvnames) % count only the names, not the interpolating functions or p or t
            if ~endsWith(mvnames{im},'_fn') & ~(strcmp(mvnames{im},'p') | strcmp(mvnames{im},'t'))
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
lFontSize = 18;
LineWidth = 2;
    if POROUS
        if ~Params.NOPLOTS
        figure(230+iT);
        if ~Params.HOLD
            clf
        end
        hold on
%         plot(Pmantle_MPa,[por_in(iT).vp]-[por_out(iT).vp])
%         plot(Pmantle_MPa,por_in(iT).vs-por_out(iT).vs,'--')    
        plot(Pmantle_MPa,por_out(iT).por*100,'LineWidth',LineWidth)
%         plot(Pmantle_MPa,(por_in(iT).den-por_out(iT).den)./por_in(iT).den*100)
%         plot(Pmantle_MPa,(por_in(iT).vp-por_out(iT).vp)./por_in(iT).vp*100)
%         plot(Pmantle_MPa,(por_in(iT).vs-por_out(iT).vs)./por_in(iT).vs*100,'--')    
        xlabel('P_{rock} (MPa)','FontSize',lFontSize)
%         ylabel('v-v_{porous} (m s^{-1})');
%         ylabel('$\frac{X-X_{porous}}{X}$ (\%)','interpreter','latex');
        ylabel('Porosity \phi (\%)','interpreter','latex','FontSize',lFontSize);
%         legend('\phi','\rho','V_P','V_S')
        box on
        axis tight
        end     
        
        if ~Params.NOPLOTS
            figure(270+iT);
        if ~Params.HOLD
            clf
        end
        hold on
%         plot(Pmantle_MPa,[por_in(iT).vp]-[por_out(iT).vp])
%         plot(Pmantle_MPa,por_in(iT).vs-por_out(iT).vs,'--')    
        hl = plot(por_out(iT).por*100,r_mantle_m*1e-3,'LineWidth',LineWidth);

%         plot((por_in(iT).den-por_out(iT).den)./por_in(iT).den*100,r_mantle_m*1e-3)
%         plot((por_in(iT).vp-por_out(iT).vp)./por_in(iT).vp*100,r_mantle_m*1e-3)
%         plot((por_in(iT).vs-por_out(iT).vs)./por_in(iT).vs*100,r_mantle_m*1e-3,'--')   
%         set(gca,'ydir','reverse');
        ylabel('r_{rock} (km)','FontSize',lFontSize)
%         ylabel('v-v_{porous} (m s^{-1})');
%         xlabel('$\frac{X-X_{porous}}{X}$ (\%)','interpreter','latex');
                xlabel('Porosity $\phi$ ($\%$)','interpreter','latex','FontSize',lFontSize);

%         legend('\phi','\rho','V_P','V_S')
        box on
        axis tight
        end
        if ~Params.NOPLOTS
            figure(280+iT);
        if ~Params.HOLD
            clf
        end
        hold on
%         plot(Pmantle_MPa,[por_in(iT).vp]-[por_out(iT).vp])
%         plot(Pmantle_MPa,por_in(iT).vs-por_out(iT).vs,'--')    
%         plot(squeeze(log10(permeability)),r_mantle_m'*ones(1,length(permeability(1,1,:)))*1e-3)
        hl1 = plot(squeeze(log10(permeability(1,:,1))),r_mantle_m*1e-3,'LineWidth',LineWidth);
        hl2 = plot(squeeze(log10(permeability(1,:,2))),r_mantle_m*1e-3,'--','LineWidth',LineWidth);
        hl3 = plot(squeeze(log10(permeability(1,:,3))),r_mantle_m*1e-3,'LineWidth',LineWidth);
        hl4 = plot(squeeze(log10(permeability(1,:,4))),r_mantle_m*1e-3,'LineWidth',LineWidth);
        hl5 = plot(squeeze(log10(permeability(1,:,5))),r_mantle_m*1e-3,'LineWidth',LineWidth);

        legend({'Crust in general','Upper crust in general','Low permeability upper crust','Disturbed crust','Oceanic crust'});

%         set(gca,'ydir','reverse');
        ylabel('r_{mantle} (km)')
%         ylabel('v-v_{porous} (m s^{-1})');
        xlabel('log$_{10}$ permeability','interpreter','latex');
        box on
        end
    end
end


%% create the table of depths and heat flux values
% multirow{4}{*}{5 Wt\%} &Ice Ih &-&141& 127  & 96 & 63  & 24 \\
% & Liquid &-&63&92&128&244&539\\
% & Ice III &-&29&-&-&-&-\\
% & Ice V &-&158&134&49&-&-\\
% & Ice VI &-&409&411&411&355&237\\

Tb_str = getTableStr(Planet.Tb_K,Planet.Tb_K); 
Qb_str = getTableStr(Planet.Tb_K,round(Qb*1e3)); 
Qc_str = getTableStr(Planet.Tb_K,round(Q_Wm2*1e3));
nu_str = getTableStr(Planet.Tb_K,log10(nu));
dzI_str = getTableStr(Planet.Tb_K,round(zI_m*1e-3));
dzOcean_str = getTableStr(Planet.Tb_K,round(dzOcean_m*1e-3));
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

if Planet.FeCore
    for iT = 1:nTbs
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
Qradrockstr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(Planet.Qmantle_Wm2*4*pi*R_Planet_m^2/1e9) '}'];
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
disp(['&$D_\mathrm{ocean}$ (km)' dzOcean_str{:} ' \\'])
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



%% Sound Speeds in the Ice and Ocean
if Params.CALC_NEW_SOUNDSPEEDS
    disp('Computing sound speeds')
    velsIce = iceVelsGagnon1990(P_MPa,T_K);
    [vfluid_kms,Ksfluid_GPa] = deal(zeros(nTbs,nsteps));
    for iT = 1:nTbs
        tic
%         out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'Ih');
%         velsIce.VIl_kms(iT,:) = 1e-3*out.Vp;
%         velsIce.VIt_kms(iT,:) = 1e-3*out.Vs;
%         velsIce.KsI_GPa(iT,:) = 1e-3*out.Ks;
%         velsIce.GsI_GPa(iT,:) = 1e-3*out.shear;
% 
%         out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'III');
%         velsIce.VIIIl_kms(iT,:) = 1e-3*out.Vp;
%         velsIce.VIIIt_kms(iT,:) = 1e-3*out.Vs;
%         velsIce.KsIII_GPa(iT,:) = 1e-3*out.Ks;
%         velsIce.GsIII_GPa(iT,:) = 1e-3*out.shear;
% 
%         out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'V');
%         velsIce.VVl_kms(iT,:) = 1e-3*out.Vp;
%         velsIce.VVt_kms(iT,:) = 1e-3*out.Vs;
%         velsIce.KsV_GPa(iT,:) = 1e-3*out.Ks;
%         velsIce.GsV_GPa(iT,:) = 1e-3*out.shear;
%         
%         out=SeaFreeze([P_MPa(iT,:)' T_K(iT,:)'],'VI');
%         velsIce.VVIl_kms(iT,:) = 1e-3*out.Vp;
%         velsIce.VVIt_kms(iT,:) = 1e-3*out.Vs;
%         velsIce.KsVI_GPa(iT,:) = 1e-3*out.Ks;
%         velsIce.GsVI_GPa(iT,:) = 1e-3*out.shear;
        
        ir = find(phase(iT,:)==0);
        vfluid_kms(iT,ir) = fluidSoundSpeeds(P_MPa(iT,ir),T_K(iT,ir),wo,Planet.Ocean.comp);
        ireal = find(~isnan(vfluid_kms(iT,:)));
        if length(ireal)<length(ir)
            disp('WARNING: extrapolating fluid sound speeds; this is seawater above 120 MPa, right?')
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
    end
    save([datpath savefile 'Vels'],'Ksfluid_GPa','velsIce','vfluid_kms');
else
    load([datpath savefile 'Vels']);
end

%% Electrical Conductivity
if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
    k_S_m = NaN*ones(size(P_MPa)); 
    switch Planet.Ocean.comp
        case 'MgSO4'
            LarionovMgSO4 = getLarionov1984MgSO4Conductivities(1);
            Pextrap = LarionovMgSO4.Pextrap_MPa;
            Textrap = LarionovMgSO4.Textrap_K;
            k_S_m_extrap = LarionovMgSO4.k_S_m_extrap_p01m;
            [Pgg,Tgg]=meshgrid(Pextrap,Textrap);
            for iTb = 1:length(Planet.Tb_K)
                k_S_m(iTb,phase(iTb,:)==0) = interp2(Pextrap,Textrap,k_S_m_extrap,P_MPa(iTb,phase(iTb,:)==0),T_K(iTb,phase(iTb,:)==0),'spline') * conduct_scaling_MgSO4;
            end
        case 'Seawater'
            if wo>2
                for iTb = 1:length(Planet.Tb_K)
                    k_S_m(iTb,phase(iTb,:)==0)=swEOS.gsw.C(wo*ones(1,length(T_K(iTb,phase(iTb,:)==0))),T_K(iTb,phase(iTb,:)==0),P_MPa(iTb,phase(iTb,:)==0)*10);
                end
            end  
    end
end

%% construct output for seismic modeling:
        % ice                           ocean                                            mantle         core 
% [g_Planet_ms2,P_Planet_MPa,T_Planet_K,r_Planet_m,rho_pPlanet_kgm3,VP_Planet_kms,VS_Planet_kms,QS_overfgamma_Planet,k_S_mPlanet]=...
%     deal(zeros(nTbs,indSil(1)-1+nsteps_mantle(1)+Params.nsteps_core));
%% Plot settings
LineWidth=1;
SoundLineWidth=1.5;
if Params.HOLD
    LineWidth = LineWidth+0.5;
    SoundLineWidth = SoundLineWidth+0.5;
end
LineStyle=Params.LineStyle;       
ymax = 1.05*R_Planet_m*1e-3;
nfig = 3000;

if Planet.FeCore
    [VP_Planet_kms,VS_Planet_kms,Ks_Planet_GPa,Gs_Planet_GPa,QS_overfgamma_Planet,k_S_mPlanet,phasePlanet] = ...
        deal(nan(nTbs,length([g_ms2(1,1:indSil(1)-1) interior(1).g_ms2 g_ms2_core(1,:)])));
else
    [VP_Planet_kms,VS_Planet_kms,Ks_Planet_GPa,Gs_Planet_GPa,QS_overfgamma_Planet,k_S_mPlanet,phasePlanet] = ...
        deal(nan(nTbs,length([g_ms2(1,1:indSil(1)-1) interior(1).g_ms2])));
end
for iT = 1:nTbs
    % grab the indices for the ices, but only for the portion of the
    % interior calculation above the silicate interface
    H2Oinds = 1:indSil(iT)-1;
    indsI = find(phase(iT,H2Oinds)==1);
    indsLiquid = find(phase(iT,H2Oinds)==0);
    indsIII = find(phase(iT,H2Oinds)==3);
    indsV = find(phase(iT,H2Oinds)==5);
    indsVI = find(phase(iT,H2Oinds)==6);
    
    phasePlanet(iT,indsI) = 1;
    phasePlanet(iT,indsLiquid) = 0;
    phasePlanet(iT,indsIII) = 3;
    phasePlanet(iT,indsV) = 5;
    phasePlanet(iT,indsVI) = 6;
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
    
    VP_Planet_kms(iT,indsI) = velsIce.VIl_kms(iT,indsI);
    VP_Planet_kms(iT,indsLiquid) = vfluid_kms(iT,indsLiquid);
    VP_Planet_kms(iT,indsIII) = velsIce.VIIIl_kms(iT,indsIII);
    VP_Planet_kms(iT,indsV) = velsIce.VVl_kms(iT,indsV);
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
    Hp_iceI = Seismic.g_aniso_iceI*T_K(iT,indsI);
    QS_overfgamma_iceI = Seismic.B_aniso_iceI*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceI./T_K(iT,indsI))/Seismic.LOW_ICE_Q;

    Tthis=T_K(iT,indsIII);
    Hp_iceIII = Seismic.g_aniso_iceIII*Tthis;
    QS_overfgamma_iceIII = Seismic.B_aniso_iceIII*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceIII./Tthis)/Seismic.LOW_ICE_Q;
    Tthis=T_K(iT,indsV);
    Hp_iceV = Seismic.g_aniso_iceV*Tthis;
    QS_overfgamma_iceV = Seismic.B_aniso_iceV*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceV./Tthis)/Seismic.LOW_ICE_Q;
    Tthis=T_K(iT,indsVI);
    Hp_iceVI = Seismic.g_aniso_iceVI*Tthis;
    QS_overfgamma_iceVI = Seismic.B_aniso_iceVI*....
        exp(Seismic.gamma_aniso_iceI*Hp_iceVI./Tthis)/Seismic.LOW_ICE_Q;

    QS_overfgamma_Planet(iT,indsI) = QS_overfgamma_iceI;
    QS_overfgamma_Planet(iT,indsLiquid) = 0*vfluid_kms(iT,indsLiquid);
    QS_overfgamma_Planet(iT,indsIII) = QS_overfgamma_iceIII;
    QS_overfgamma_Planet(iT,indsV) = QS_overfgamma_iceV;
    QS_overfgamma_Planet(iT,indsVI) = QS_overfgamma_iceVI;
    QS_overfgamma_Planet(iT,indSil(iT):indSil(iT)+length(interior(iT).VS_mantle_kms)-1) = interior(iT).QS_overfgamma;
    if Planet.FeCore
        QS_overfgamma_Planet(iT,start_core:start_core-1+length(thisVPcore)) = thisQScore;
    end

 %%save the data to a text file
    Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [P_Planet_MPa(iT,:)', T_Planet_K(iT,:)', r_Planet_m(iT,:)'*1e-3, rho_pPlanet_kgm3(iT,:)',VP_Planet_kms(iT,:)',VS_Planet_kms(iT,:)',QS_overfgamma_Planet(iT,:)' Ks_Planet_GPa(iT,:)' Gs_Planet_GPa(iT,:)' g_Planet_ms2(iT,:)' phasePlanet(iT,:)'];
        header = sprintf('%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',...
            'P (MPa)','T (K)','r (km)','rho (kg m-3)','VP (km s-1)','VS (km s-1)','QS/gamma','KS (GPa)','GS (GPa)','g (m/s2)','phase');
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY %11/8/19 added g and phase. This changes the column number of electrical conductivity
%         k_S_mPlanet(iT,:) = [0*QS_overfgamma_iceI(iT,:) k_S_m(iT,indsLiquid)...
%             0*QS_overfgamma_iceIII(iT,indsIII) 0*QS_overfgamma_iceV(iT,indsV) 0*QS_overfgamma_iceVI(iT,indsVI) ...
%             0*interior(iT).QS_overfgamma thiskScore];
        k_S_mPlanet(iT,indsI) = 0*QS_overfgamma_iceI;
        k_S_mPlanet(iT,indsLiquid) = k_S_m(iT,indsLiquid);
        k_S_mPlanet(iT,indsIII) = 0*QS_overfgamma_iceIII;
        k_S_mPlanet(iT,indsV) = 0*QS_overfgamma_iceV;
        k_S_mPlanet(iT,indsVI) = 0*QS_overfgamma_iceVI;
        k_S_mPlanet(iT,indSil(iT):indSil(iT)+length(interior(iT).VS_mantle_kms)-1) = 0*interior(iT).QS_overfgamma;
        if Planet.FeCore
            k_S_mPlanet(iT,start_core:start_core-1+length(thisVPcore)) = 0*thisQScore;
        end
        Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma k_S_mPlanet(iT,:)' ]; %#ok<AGROW>
        header = sprintf('%s%s\t%s\t',header,'k (S m-1)');
    end
    thissavestr = [savefile 'Ts' num2str(Planet.Tsurf_K,'%0.0f') 'Zb' strLow num2str(Zb2(iT),'%0.0f') ...
        'mQm' num2str(1000*Planet.Qmantle_Wm2(iT),'%0.0f') 'mWm2_CMR2p' ...
        num2str(10000*Planet.Cmeasured,'%0.0f') '_' thiseos];
%     save(thissavestr,'Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma','-ascii');
    dlmwrite([datpath thissavestr '.txt'],header,'delimiter','');
    dlmwrite([datpath thissavestr '.txt'],Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma,...
             'delimiter','\t',...
             'precision','%3.5e',...
             '-append');

         dlmwrite([datpath thissavestr '_mantle.dat'],[10*interior(iT).Pmantle_MPa' interior(iT).Tmantle_K'],'delimiter','\t');


%% plot the seismic data and attenuation
if ~Params.NOPLOTS
    figure(nfig+iT);
    if ~Params.HOLD
        clf
    end
    set(gcf,'Position', [476   662   560   220],'Name',thissavestr)
    hp = subplot(1,3,1);
    if Params.HOLD
        hold on
    end
    plot(VS_Planet_kms(iT,:)',r_Planet_m(iT,:)'*1e-3,...
        VP_Planet_kms(iT,:)',r_Planet_m(iT,:)'*1e-3,'--','LineWidth',LineWidth)
        set(gcf,'color','w')
        xlabel('Sound Speeds (km s^{-1})')
        ylabel(['r_{' Planet.name '} (km)']);
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
        'LineWidth',LineWidth);
        xlabel('P (bar), T (K), \rho (kg m^{-3})');
        xmax = max([max(rho_pPlanet_kgm3(iT,:)) max(10*P_Planet_MPa(iT,:)) max(T_Planet_K(iT,:))]);
    else
        plot(P_Planet_MPa(iT,:)',r_Planet_m(iT,:)'*1e-3,...
            T_Planet_K(iT,:)',r_Planet_m(iT,:)'*1e-3,'--',...
            rho_pPlanet_kgm3(iT,:)',r_Planet_m(iT,:)'*1e-3,'k-.',...
        'LineWidth',LineWidth);
        xlabel('P (MPa), T (K), \rho (kg m^{-3})');
        xmax = max([max(rho_pPlanet_kgm3(iT,:)) max(P_Planet_MPa(iT,:)) max(T_Planet_K(iT,:))]);
    end
    set(gca,'ylim',[0 ymax],'xlim',[0 xmax]);
    grid on
    
    subplot(1,3,3)
    if Params.HOLD
        hold on
    end
    hp = plot(QS_overfgamma_Planet(iT,:)',r_Planet_m(iT,:)'*1e-3,'LineWidth',LineWidth);
    set(gca,'xscale','log','ylim',[0 ymax],'xlim',[10 1e7],'XTick',[10 100 1000 1e4 1e5 1e6 1e7])
    grid on
    xlabel('Q_S/\omega^{\gamma}')
    set(gcf,'color','w')
    
    try
        saveas(gcf,[figpath thissavestr 'QS'],Params.savefigformat);
    catch
        error(['Couldn''t save file ' thissavestr ' to the ''figures'' directory. Try setting active folder to PlanetProfile directory in Matlab.']);
    end
    if ~Params.NOPLOTS
    figure(nfig+iT+50);
    if ~Params.HOLD
        clf;
    end
    hold on
    hp = plot(Ks_Planet_GPa(iT,:)',r_Planet_m(iT,:)'*1e-3,...
        Gs_Planet_GPa(iT,:)',r_Planet_m(iT,:)'*1e-3,'--','LineWidth',LineWidth);
    set(hp,'LineWidth',LineWidth);
    ylabel(['r_{' Planet.name '} (km)'],'FontSize',lFontSize);
    xlabel('G_S and K_S (GPa)','FontSize',lFontSize);
    axis tight
    try
        saveas(gcf,[figpath thissavestr 'GsKs'],Params.savefigformat);
    catch
        error(['Couldn''t save file ' thissavestr ' to the ''figures'' directory. Try setting active folder to PlanetProfile directory in Matlab.']);
    end
   end
end
end

if ~Params.NOPLOTS

    revcop =    colormap('copper');
revcop = revcop(256:-1:1,:);
colormap(revcop);
inferno_data = CM_inferno;
colormap(inferno_data);

figure(1111);clf; clear opts
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
subplot(2,3,1); plotSolidInterior('rho','Density (kg m^{-3})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
subplot(2,3,2); plotSolidInterior('cp','Cp (J m^{-3} K^{-1})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
subplot(2,3,3); plotSolidInterior('alpha','\alpha (K^{-1})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
subplot(2,3,4); plotSolidInterior('vp','V_P (km s^{-1})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
subplot(2,3,5); plotSolidInterior('vs','V_S (km s^{-1})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)

if isfield(mantle,'Ks')
    subplot(2,3,6); plotSolidInterior('Ks','K_S (GPa)',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
    
    figure(1112); clf;
    opts.Punits = 'GPa';
    colormap(inferno_data);
        subplot(2,2,1); plotSolidInterior('rho','Density (kg m^{-3})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)        
        subplot(2,2,2); plotSolidInterior('cp','Cp (J m^{-3} K^{-1})',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
        subplot(2,2,3); plotSolidInterior('Ks','K_S (GPa)',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
        subplot(2,2,4); plotSolidInterior('Gs','G_S (GPa)',T_Planet_K,P_Planet_MPa,mantle,core,tcore,pcore,opts)
end
end

% save the mineral compositions plot 
if ~Params.NOPLOTS
figure(1111)
saveas(gcf,[figpath thissavestr],Params.savefigformat);
end

%% create the legend that describes tb, z_b, and heat flux
lstr_2 = {};
for iT = 1:nTbs
    lstr_2 = {lstr_2{:} ['T_{b}:' num2str(Planet.Tb_K(iT),'%0.2f') ' K, q_{b}/q_{c}:'...
        num2str(1e3*Qb(iT),'%0.0f') '/' num2str(1e3*Q_Wm2(iT),'%0.0f') ' mW m^{-2}, z_{b}:' num2str(1e-3*Zb2(iT),'%0.0f') ' km']};
end

figure(1122);clf;
for iT = 1:nTbs
    subplot(1,nTbs,iT); hold on
    input.phase = phasePlanet(iT,:);
    input.r_m = r_Planet_m(iT,:);
    PlanetWedge(input);
    box off
    if iT==1
        ylabel(Planet.name,'FontSize',24)
    else
      set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none')
    end
    title(['T_b = ' num2str(Tb_K(iT)) ' K'],'FontSize',16,'FontWeight','normal')
end
    

%%  plot profile with 4 subplots
%% 
if ~Params.NOPLOTS
if Params.foursubplots
Dsil_km=(R_Planet_m-R_sil_mean_m(:))*1e-3;

figure(229);
maxScale = 1.01;
minScale = 0.99;
aFontSize = 16;
lFontSize = 18;
if ~Params.HOLD
    clf;
end
%     set(gcf,'Position', [411    52   764   468])
if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
    subplot(3,6,4:6);
else
    subplot(2,6,4:6);
end
    hold on
for iT = 1:nTbs
    line(T_K(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),...
        'LineWidth',LineWidth,'LineStyle',LineStyle);
    hm = line(T_K(iT,indSil(iT)),z_m(iT,indSil(iT))*1e-3,'Color',Params.colororder(iT),'Marker','o');
    if strcmp(LineStyle,'-')
        hm.MarkerFaceColor=Params.colororder(iT);
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
else
    set(gca,'ydir','reverse',...
        'xlim',[245 max(max(T_K))],'ylim',[0 maxScale*max(Dsil_km)],...
        'FontSize',aFontSize,'YAxisLocation','right');%'XAxisLocation','top');
end
xlabel('Temperature (K)','FontSize',lFontSize);
ylabel('Depth (km)','FontSize',lFontSize);
% for iT=1:nTbs
%     hline(Dsil_km(iT),Params.colororder(iT));
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
    hp(iT) =  line(vfluid_kms(iT,Params.nsteps_iceI+1:indSil(iT)),z_m(iT,Params.nsteps_iceI+1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
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
   ax.FontSize = aFontSize;

   box on

   if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
    subplot(3,6,11); 
else
    subplot(2,6,11);
end
hold on;
for iT = 1:nTbs
    line(velsIce.VIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
       maxV(iT) = max([velsIce.VIt_kms(iT,1:indSil(iT)) ...
            velsIce.VIIt_kms(iT,1:indSil(iT)) ...
            velsIce.VIIIt_kms(iT,1:indSil(iT)) ...
            velsIce.VVt_kms(iT,1:indSil(iT)) ...
            velsIce.VVIt_kms(iT,1:indSil(iT))]);
       minV(iT) = min([velsIce.VIt_kms(iT,1:indSil(iT)) ...
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
   ax.FontSize = aFontSize;
   ax.YTickLabel = [];
box on
xlabel('Sound Speed (km s^{-1})','FontSize',lFontSize);


if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
    subplot(3,6,12); 
else
    subplot(2,6,12);
end
hold on;
for iT = 1:nTbs
    line(velsIce.VIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    maxV(iT) = max([velsIce.VIl_kms(iT,1:indSil(iT)) ...
            velsIce.VIIl_kms(iT,1:indSil(iT)) ...
            velsIce.VIIIl_kms(iT,1:indSil(iT)) ...
            velsIce.VVl_kms(iT,1:indSil(iT)) ...
            velsIce.VVIl_kms(iT,1:indSil(iT))]);
    minV(iT) = min([velsIce.VIl_kms(iT,1:indSil(iT)) ...
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
   ax.FontSize = aFontSize;
      ax.YAxisLocation = 'right';


   ylabel('Depth (km)','FontSize',lFontSize);


box on
end
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
     if ~Params.NOPLOTS  
        if wo>0
        subplot(3,6,16:18);hold on
        switch Planet.Ocean.comp
            case 'MgSO4'
                for iT = 1:nTbs                 
                    line(k_S_m(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,...
                        'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',LineWidth);
                    mink_val(iT) = min(k_S_m(iT,1:indSil(iT)));
                    maxk_val(iT) = max(k_S_m(iT,1:indSil(iT)));
                end
            case 'Seawater'
                for iT = 1:nTbs
                    line(k_S_m(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),...
                        'LineStyle',LineStyle,'LineWidth',LineWidth);
                    mink_val(iT) = min(k_S_m(iT,1:indSil(iT)));
                    maxk_val(iT) = max(k_S_m(iT,1:indSil(iT)));
                end
            case 'NH3'
                for iT = 1:nTbs
%                     line(k_S_m(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',LineWidth);
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
        ax.FontSize = aFontSize;
        ax.YAxisLocation = 'right';

        xlabel('Electrical Conductivity (S m^{-1})','FontSize',aFontSize);
        ylabel('Depth (km)','FontSize',aFontSize);
        box on
        end
    end
end


% subplot(3,2,[1 3 5]);
if ~Params.NOPLOTS
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
    ht(iT)= plot(P_MPa(iT,1:indSil(iT)),rho_kgm3(iT,1:indSil(iT)),Params.colororder(iT),...
        'LineWidth',LineWidth,'LineStyle',LineStyle); 
    hm = plot(Psil_MPa(iT),interp1(P_MPa(iT,:),rho_kgm3(iT,:),Psil_MPa(iT)),[Params.colororder(iT) 'o']);
    if ~Params.HOLD
        hm.MarkerFaceColor=Params.colororder(iT);
    end
end

for ir = 1:length(Params.wref)
    hw(ir) = plot(Pref_MPa,rho_ref_kgm3(ir,:),['k' Params.wrefLine]);
end

if Params.Legend
    hleg1 = legend([ht],lstr_2{:});
    set(hleg1,'location',Params.LegendPosition,'box','off')
end

box on;
ax = gca;
ax.FontSize = aFontSize;
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
xlabel('Pressure (MPa)','FontSize',lFontSize)
ylabel('Density (kg m^{-3})','FontSize',lFontSize)

    try
        saveas(gcf,[figpath savefile],Params.savefigformat);
    catch
        error(['Couldn''t save file ' savefile ' to the ''figures'' directory. Maybe it doesn''t exist in the PP folder?']);
    end
end

%% plot profile with 2 subplots
else
if ~Params.NOPLOTS
figure(229);
maxScale = 1.01;
if ~Params.HOLD
    clf;
end
    set(gcf,'Position', [411    52   764   468])
    subplot(1,2,2);
    hold on

Dsil_km=(R_Planet_m-R_sil_mean_m(:))*1e-3;

for iT = 1:nTbs
    hp(iT) =  line(vfluid_kms(iT,Params.nsteps_iceI+1:indSil(iT)),z_m(iT,Params.nsteps_iceI+1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VIIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VIIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VVl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VVt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VVIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
    line(velsIce.VVIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','-.','LineWidth',SoundLineWidth);
end
for iT = 1:nTbs
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        line(k_S_mMgSO4p01Planet(iT,1:indSil(iT))*25,z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle','--','LineWidth',LineWidth);
    end
end
set(gca,'ydir','reverse','xlim',[0 5],'ylim',[0 maxScale*max(Dsil_km)]);%,'xlim',[1 4]);
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        xlabel('Conductivity (S m^{-1} \times 25) and Sound Speed (km s^{-1})');
    else
        xlabel('Sound Speed (km s^{-1})');
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
    line(ax2,T_K(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineWidth',LineWidth);
    hm = line(T_K(iT,indSil(iT)),z_m(iT,indSil(iT))*1e-3,'Color',Params.colororder(iT),'Marker','o');
    hm.MarkerFaceColor=Params.colororder(iT);
end

set(gca,'ydir','reverse','xlim',[245 320],'ylim',[0 maxScale*max(Dsil_km)]);
xlabel('Temperature (K)');
ylabel('Depth (km)');
% for iT=1:nTbs
%     hline(Dsil_km(iT),Params.colororder(iT));
% end
box on;



subplot(1,2,1);
hold all 
R2ind = C2mean+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
Psil_MPa = diag(P_MPa(:,R2ind(:)));
for iT=1:nTbs
    ht(iT)=  plot(P_MPa(iT,1:indSil(iT)),rho_kgm3(iT,1:indSil(iT)),Params.colororder(iT),'LineWidth',LineWidth); 
    hm = plot(Psil_MPa(iT),interp1(P_MPa(iT,:),rho_kgm3(iT,:),Psil_MPa(iT)),[Params.colororder(iT) 'o']);
    hm.MarkerFaceColor=Params.colororder(iT);
end

hw(1) = plot(Pref_MPa,rho_ref_kgm3(1,:),'k--');
hw(2) = plot(Pref_MPa,rho_ref_kgm3(2,:),'k--');
hw(3) = plot(Pref_MPa,rho_ref_kgm3(3,:),'k--');
hw(4) = plot(Pref_MPa,rho_ref_kgm3(4,:),'k--');


% hleg1 = legend([ht hw],{lstr_2{:},'0 wt%','5 wt%','10 wt%','15 wt%'});%,'20 wt%');
if Params.Legend
    hleg1 = legend([ht],lstr_2{:});
    set(hleg1,'location','southeast','box','off')
end

axis tight; box on;
set(gca,'xlim',[0 maxScale*max(Psil_MPa)])
xlabel('Pressure (MPa)')
ylabel('Density (kg m^{-3})')

% ax(1) = gca;
% ax(2)=axes('Position',get(ax(1),'Position'),...
%    'XAxisLocation','top',...
%    'YAxisLocation','right',...
%    'XColor','none',...
%    'Color','none');
% 
% % add the sound speeds
% for iT = 1:nTbs
%     line(P_MPa(iT,:),vfluid_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VIl_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VIt_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VIIl_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VIIt_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VIIIl_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VIIIt_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VVl_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VVt_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VVIl_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
%     line(P_MPa(iT,:),velsIce.VVIt_kms(iT,:),'Color',Params.colororder(iT),'Parent',ax(2),'LineStyle','-.','LineWidth',SoundLineWidth);
% end
% set(ax(2),'ylim',[1.5 5])
% set(ax(2),'xlim',[0 1800])
% ylabel('Sound Speeds (km s^{-1})')

    saveas(gcf,[figpath savefile],Params.savefigformat);
end
end
%%
if ~Params.NOPLOTS
figure(228)
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
% plot(g_ms2',r_m'*1e-3,gp,R_sil_mean_m*1e-3,'o');
hl = plot(g_Planet_ms2',r_Planet_m'*1e-3);
hold on
for iT = 1:nTbs
    hp = plot(gp(iT),R_sil_mean_m*1e-3,'o');
    set(hl(iT),'LineWidth',LineWidth,'LineStyle',LineStyle,'Color',Params.colororder(iT));
    set(hp,'Color',Params.colororder(iT));
end
xlabel('g (m s^{-2})')
ylabel('R (km)')
axis tight

subplot(1,2,2);
if Params.HOLD
    hold on;
end
hl = plot(P_Planet_MPa',r_Planet_m'*1e-3);
for iT = 1:nTbs
    set(hl(iT),'LineWidth',LineWidth,'LineStyle',LineStyle,'Color',Params.colororder(iT));
end

 xlabel('Pressure (MPa)')
 ylabel('r_{Planet} (km)')
axis tight

%%    
% figure(232);clf;
% subplot(2,1,2)
% hold all;    
% clear ht hw
% for iT=1:nTbs
%     ht(iT)=  plot(P_MPa(iT,:),T_K(iT,:),Params.colororder(iT));
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
end
end %PlanetProfile
%% properties of water ice
%Note: for VspChoukroun inds: liquid = 1, Ice I = 2, II = 3, III = 4, V = 5, VI = 6
function rho_kgm3 = getRhoIce(P_MPa,T_K,ind)
%     % convert to appropriate phase using our nomenclature
    if ~(ind==5 || ind==6)
        ind = ind+1;
    end
    switch ind
        % 1 (liquid water) isn't possible because of the above convention used for
        % implementing Choukroun and Grasset's (2010) EOS.
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
    end
    try
        out = SeaFreeze([P_MPa,T_K],material);
        rho_kgm3 = out.rho;
    catch
        if ind<4 % ices I and II
            rho_kgm3 = 1000./getVspChoukroun2010(P_MPa,T_K,ind);
        else
            MW_H2O = 18.01528;
            switch ind
                case 4 % ice III
                    iceIII_Vfn = load('iceIII_sp_V_fPT');
                    rho_kgm3 = 1./fnval(iceIII_Vfn.sp_V_fPT,{1e-3*P_MPa T_K});
                case 5 % ice V
                    iceV_Vfn = load('iceV_sp_V_fPT');
                    rho_kgm3 = 1./fnval(iceV_Vfn.sp_V_fPT,{1e-3*P_MPa T_K});            
                case 6 % ice VI
                    rho_kgm3 = MW_H2O*1000./iceVI_PT_EOS_bezacier(1e-3*P_MPa,T_K);
            end
        end
    end
end %getRhoIce
function [Cp,alpha] = getCpIce(P_MPa,T_K,ind)
    % convert to appropriate phase using our nomenclature
    if ~(ind==5 || ind==6)
        ind = ind+1;
    end
    switch ind
        % 1 (liquid water) isn't possible because of the above convention used for
        % implementing Choukroun and Grasset's (2010) EOS.
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
    end
    try
        out = SeaFreeze([P_MPa,T_K],material);
        Cp = out.Cp;
        alpha = out.alpha;
    catch
        disp(['T_ice = ' num2str(T_K) '. This seems to be too low for SeaFreeze. Using Choukroun and Grasset (2010) instead.']);
        Cp = CpH2O_Choukroun(P_MPa,T_K,ind);
        alpha = [];
    end
end %getCpIce     
function phase = getIcePhase(P_MPa,T_K,w_pct,str_comp)
    switch str_comp
        case 'MgSO4' 
%             if P_MPa<800
%                 W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
%                 mo=1000./W_MgSO4./(1./(0.01.*w_pct)-1); %conversion to molality from Wt%
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
            phase = getIcePhaseNH3(P_MPa,T_K,w_pct);
    end
end % getIcePhase
function Tfreeze_K = getTfreeze(P_MPa,wo,str_comp,Tprior)
    switch str_comp
        case 'MgSO4'
%             if P_MPa<800
%                 W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
%                 mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from Wt%
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
                Tfreeze_K = fzero(@(T) 0.5-LBFIcePhase(P_MPa,T,wo,'NaCl'),Trange,options);
        case 'Seawater'
            global swEOS
            Tfreeze_K = swEOS.gsw.tfreezing(wo,10*P_MPa);
        case 'NH3'
            Tfreeze_K = fzero(@(T) L_IceNH3(P_MPa,T,wo,1),[200 350]);
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
            Pfreeze_MPa = fzero(@(P) 0.5-LBFIcePhase(P,T_K,wo,'NaCl'),[0.1 209],options); % the P bounds are sensitive because multivalued solutions (ices Ih and III) will cause fzero to fail
        case 'Seawater'
            global swEOS
            Pfreeze_MPa = 0.1*swEOS.gsw.pfreezing(wo,T_K);
        case 'NH3'
            Pfreeze_MPa = fzero(@(P) L_IceNH3(P,T_K,wo,1),[0 250]);
    end
end %getPfreeze
function Pfreeze_MPa = getPfreezeIII(T_K,wo,str_comp)
    switch str_comp
       case 'MgSO4'             
           Pfreeze_MPa = fzero(@(P) L_IceMgSO4(P,T_K,wo,3),[0 500]);
       case 'NaCl'
%            Pfreeze_MPa = fzero(@(P) 0.5-LBFIcePhase(P,T_K,wo,'NaCl'),[0.1 500]);
       case 'NH3'
           Pfreeze_MPa = fzero(@(P) L_IceNH3(P,T_K,wo,3),[0 500]);
    end
end %getPfreezeIII
%% fluid properties
function [rho_kgm3,Cp,alpha_Km1]=fluidEOS(P_MPa,T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'
            W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
            mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from Wt%
            [rho_kgm3,Cp,alpha_Km1]=MgSO4_EOS2_planetary_smaller(mo,P_MPa,T_K-273.15);
        case 'NaCl'
            global swEOS
            disp('WARNING: NaCl is not yet implemented for fluidEOS.')
        case 'Seawater'
            global swEOS extrapflag
            % expect that wo is absolute salinity
            rho_kgm3=swEOS.gsw.dens(wo,T_K,P_MPa*10);
            Cp = swEOS.gsw.cp(wo,T_K,P_MPa*10);
            alpha_Km1 = swEOS.gsw.alpha(wo,T_K,P_MPa*10);
            if isnan(rho_kgm3)
                if ~extrapflag % only warn once
                    disp('WARNING: extrapolating fluid rho, Cp, and alpha; this is seawater above 120 MPa, right?')
                    extrapflag = 1;
                end
                Pin = [110 115 120];
                rhoin=swEOS.gsw.dens(ones(1,3)*wo,ones(1,3)*T_K,Pin*10);
                Cpin = swEOS.gsw.cp(ones(1,3)*wo,ones(1,3)*T_K,Pin*10);
                alphain = swEOS.gsw.alpha(ones(1,3)*wo,ones(1,3)*T_K,Pin*10);
                
                rho_kgm3 = interp1(Pin,rhoin,P_MPa,'linear','extrap');
                Cp = interp1(Pin,Cpin,P_MPa,'linear','extrap');
                alpha_Km1 = interp1(Pin,alphain,P_MPa,'linear','extrap');                
            end
        case 'NH3'
           disp('WARNING: NH3 is not yet implemented for fluidEOS.')
    end
end %fluidEOS
function [vel_kms,Ks_GPa] = fluidSoundSpeeds(P_MPa,T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'
            W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
            mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from Wt%
            vel_kms = MgSO4_EOS2_planetary_velocity_smaller_vector(mo,P_MPa,T_K-273.15); % assumes P and T are the same length.
        case 'NaCl'
            global swEOS
            disp('WARNING: NaCl is not yet implemented for fluidSoundSpeeds.')
        case 'Seawater'
            global swEOS
            % expect that wo is absolute salinity
            vel_kms=1e-3*swEOS.gsw.vel(wo*ones(1,length(T_K)),T_K,P_MPa*10);
        case 'NH3'
            disp('WARNING: NH3 is not yet implemented for fluidSoundSpeeds.')
    end
end % fluidSoundSpeeds
function zero_alpha = alphaAdjust(P_MPa,T_K,wo,comp)
    [~,~,zero_alpha]= fluidEOS(P_MPa,T_K,wo,comp); 
end % alphaAdjust
%% Core Size
function [C2MR2,R_fe] = CoreSize(rho_s,rho_fe,C_H2O,M_above_kg,R_s)
    global M_Planet_kg R_Planet_m
    try
        R_fe = fzero(@(R_fe) getR_fe(rho_s,rho_fe,M_above_kg,R_s,R_fe),[0 R_Planet_m]);
    catch
        R_fe = NaN;
    end
    C2MR2 = C_H2O+8/15*pi*((R_s.^5-R_fe.^5)*rho_s+rho_fe*R_fe.^5);%/M_Planet_kg/R_Planet_m^2;
end %CoreSize
function zero_me = getR_fe(rho_s,rho_fe,M_above_kg,R_s,R_fe)
    global M_Planet_kg 
    zero_me = 4*pi/3*rho_fe*R_fe.^3 - (M_Planet_kg - M_above_kg-4*pi/3*rho_s*(R_s.^3-R_fe.^3));
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

if strcmp(prop,'Ks') || strcmp(prop,'Gs')
    pcolor(mantle.t,mantle.p*1e3*xp,mantle.(prop)*1e-4);
else
    pcolor(mantle.t,mantle.p*1e3*xp,mantle.(prop));
end

hold on; 
if ~isempty(core)
    propcore = core.([prop '_fn'])(pcore,tcore);
    if strcmp(prop,'Ks') || strcmp(prop,'Gs')
        pcolor(tcore,pcore*xp,propcore*1e-4);   
    else
        pcolor(tcore,pcore*xp,propcore);
    end
end
hp = plot(T_Planet_K',P_Planet_MPa'*xp,'w-','LineWidth',1);
shading interp; colorbar; 
box on; set(gca,'ydir','reverse'); 
if opts.Ttight
    xlims = get(gca,'XLim');
    set(gca,'XLim',[xlims(1) max(max(T_Planet_K))*1.01])
end
ylabel(['Pressure (' opts.Punits ')']);xlabel('Temperature (K)');
title(title_str)

end %plotSolidInterior