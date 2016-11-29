function PlanetProfile(Planet,Seismic,Params)
%% PlanetProfile
%
% S. Vance, M. Bouffard, M. Choukroun, and C. Sotin. 
% Ganymede's internal structure including thermodynamics of magnesium sulfate oceans in contact with ice. 
% Planetary And Space Science, 96:62?70, 2014.    
% (http://dx.doi.org/10.1016/j.pss.2014.03.011)


%% globals for functions in fzero
global M_Planet_kg R_Planet_m
M_Planet_kg = Planet.M_kg;
R_Planet_m = Planet.R_m;

savefile = [Planet.name 'Profile_' Planet.Ocean.comp ...
    num2str(Planet.Ocean.w_ocean_pct) 'WtPct'];

Gg = 6.67300e-11; % m3 kg-1 s-2

n_iceI=Params.nsteps_iceI;
n_ocean = Params.nsteps_ocean;

%%
if Params.CALC_NEW
nTbs = length(Planet.Tb_K);

nsteps = n_iceI+n_ocean;
[T_K,P_MPa,rho_kgm3] = deal(zeros(nTbs,nsteps));

phase = zeros(nTbs,nsteps);
phase(:,1:n_iceI)=1;
wo = Planet.Ocean.w_ocean_pct;
Tb_K = Planet.Tb_K;

%%
%ICE Ih
%--------------------------------------------------------------------------
for kt = 1:nTbs %draw four thermal profiles corresponding to four different choices of temperature at the bottom of the Ice I shell
    T_K(kt,1) = Planet.Tsurf_K;
    P_MPa(kt,1) = Planet.Psurf_MPa;

    rho_kgm3(kt,1) = getRhoIce(P_MPa(kt,1),T_K(kt,1),1);
    Pb_MPa(kt) = getPfreeze(Tb_K(kt),wo,Planet.Ocean.comp);

    deltaP = Pb_MPa(kt)/(n_iceI-1);

for il=2:n_iceI % propagate P,T,rho to the bottom of the ice
    P_MPa(kt,il) = P_MPa(kt,il-1) + deltaP;   
    T_K(kt,il) = (Planet.Tb_K(kt).^(P_MPa(kt,il)./Pb_MPa(kt))).*(Planet.Tsurf_K.^(1-P_MPa(kt,il)./Pb_MPa(kt)));
    rho_kgm3(kt,il) = getRhoIce(P_MPa(kt,il),T_K(kt,il),1); 
end
   
%OCEAN + ICE III/V/VI SHELL
%--------------------------------------------------------------------------
  deltaP = (Params.Pseafloor_MPa-Pb_MPa(kt))/n_ocean; %
for il =1:n_ocean      
    ill = il+n_iceI;
    disp(['kt: ' num2str(kt) '; il: ' num2str(il) '; P_MPa: ' num2str(round(P_MPa(kt,ill-1)))]);
    P_MPa(kt,ill) = Pb_MPa(kt) + il*deltaP;

    if il ==1 % establish the phase vector
        phase(kt,ill) = getIcePhase(P_MPa(kt,ill),T_K(kt,ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
    elseif phase(kt,ill-1)~=6
        phase(kt,ill) = getIcePhase(P_MPa(kt,ill),T_K(kt,ill-1),wo,Planet.Ocean.comp);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
        if phase(kt,ill-1)==3 && phase(kt,ill)==0 % fix to instabilities in the phase lookup for ice III
            phase(kt,ill)=3;
        end
    else
        phase(kt,ill) = 6;
    end
    
    if phase(kt,ill) == 0 %if ocean
      [rho_ocean,Cp,alpha_o]= fluidEOS(P_MPa(kt,ill),T_K(kt,ill-1),wo,Planet.Ocean.comp);
      T_K(kt,ill) = T_K(kt,ill-1)+ alpha_o*T_K(kt,ill-1)./(Cp)./(rho_ocean)*deltaP*1e6; %adiabatic gradient in ocean; this introduces an error, in that we are using the temperature from the previous step
      rho_kgm3(kt,ill) = fluidEOS(P_MPa(kt,ill),T_K(kt,ill),wo,Planet.Ocean.comp);    
    else %if ice
        if ~isfield(Planet.Ocean,'fnTfreeze_K')
            T_K(kt,ill) = getTfreeze(P_MPa(kt,ill),wo,Planet.Ocean.comp); %find the temperature on the liquidus line corresponding to P(k,i+nIceI); should really use a conductive profile here, but that would seem to naturally bring us back to the liquidus. Steve wants to confirm.
        else
            T_K(kt,ill) = Planet.Ocean.fnTfreeze_K(P_MPa(kt,ill),wo);
        end
      rho_kgm3(kt,ill) = getRhoIce(P_MPa(kt,ill),T_K(kt,ill),phase(kt,ill));
    end
end
rho_kgm3(kt,1) = rho_kgm3(kt,2); %continuity
    save(savefile,'P_MPa','Pb_MPa','T_K','Tb_K','phase','deltaP','wo','nTbs','rho_kgm3','rho_ocean','Cp','alpha_o','nsteps'); % save the progress at each step
end  
else
    load(savefile);
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
           parfor il=1:nPr
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
        end
    save([str_ref_densities],'rho_ref_kgm3','Pref_MPa','Tref_K')
else
    load([str_ref_densities]);
end

%%%%%%%%%%%%%%%%%%%%%
% convert to depth
%%%%%%%%%%%%%%%%%%%
%% calculate gravity in each layer instead of assuming surface gravity applies.
% allocate variables
    z_m = zeros(nTbs,nsteps);
    g_ms2 = z_m;
    [M_above_kg,M_below_kg] = deal(rho_kgm3);
    [z_m,r_m] = deal(zeros(nTbs,nsteps)); 
    r_m(:,1) = Planet.R_m;
    g_ms2(1:nTbs,1) = Planet.gsurf_ms2; 

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
D_conductivity = 632; % W m-1; Andersson et al. 2005 (For comparison, Mckinnon 2006 uses a value of 621 from Slack 1980)
nTbs = length(T_K(:,1));
for kt = 1:nTbs
    Qb(kt) = D_conductivity*log(Planet.Tb_K(kt)/Planet.Tsurf_K)/Zb2(kt);
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

    % figure out the indices for different ice layers
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

dindsIII = zIII_m>0;
dzOcean_m(dindsIII) = zIII_m(dindsIII)-zI_m(dindsIII);

dindsV = zV_m>0;
dzOcean_m(~dindsIII & dindsV) = zV_m(~dindsIII & dindsV)-zI_m(~dindsIII & dindsV); 
dzIII_m(dindsIII & dindsV) = zV_m(dindsV & dindsIII)- zIII_m(dindsV & dindsIII);
dindsVI = zVI_m>0;
dzOcean_m(dindsVI & ~dindsV) = zVI_m(dindsVI & ~dindsV) - zI_m(dindsVI & ~dindsV);
dzV_m(dindsV) = zVI_m(dindsV) - zV_m(dindsV);

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
        C2mean(iT) = round(mean(C2inds{iT}));
        C2max(iT) = max(C2inds{iT});
        C2min(iT) = min(C2inds{iT});
        R_sil_mean_m(iT) = R_sil_m(iT,C2mean(iT));
    end
    R_Fe_mean_m = zeros(1,nTbs);
   
    dzOcean_m(~dindsVI & ~dindsV) = R_Planet_m - R_sil_mean_m(~dindsVI & ~dindsV)-zI_m(~dindsVI & ~dindsV);

    
   % Plot the results
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
else % With a core
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
    figure(2296);clf;hold all
    for iT=1:nTbs
        C2inds{iT} = find(C2(iT,:)/MR2>Planet.Cmeasured-Planet.Cuncertainty & C2(iT,:)/MR2<Planet.Cmeasured+Planet.Cuncertainty);
        C2mean(iT) = round(mean(C2inds{iT}));
        C2max(iT) = max(C2inds{iT});
        C2min(iT) = min(C2inds{iT});
        R_Fe_mean_m(iT) = R_Fe_m(iT,C2mean(iT));
        R_sil_mean_m(iT) = R_sil_m(iT,C2mean(iT));
        R_Fe_range_m(iT) = R_Fe_m(iT,C2max(iT))-R_Fe_m(iT,C2min(iT));
        R_sil_range_m(iT) = R_sil_m(iT,C2min(iT))-R_sil_m(iT,C2max(iT));
        plot(R_Fe_m(iT,C2inds{iT})'*1e-3,R_sil_m(iT,C2inds{iT})'*1e-3);
    end

    lstr_3 = {};
    for iT = 1:nTbs
        lstr_3 = {lstr_3{:} ['T_{b}:' num2str(Planet.Tb_K(iT),'%0.1f') 'K']};
    end
    legend(lstr_3)
    box on
    xlabel('R_{Fe} (km)');
    ylabel('R_{sil} (km)')
    title(['Fe core ; C/MR2=' num2str(Planet.Cmeasured) '\pm' num2str(Planet.Cuncertainty) '; W =' num2str(wo) ' Wt%; \rho_{sil}: ' num2str(Planet.rho_sil_withcore_kgm3,'%0.0f') '; \rho_{Fe}: ' num2str(rho_Fe_kgm3,'%0.0f')])

    % Display Results For the Calculation of Mantle and Core Depths
    % use R_sil_mean_m to calculate thickness of ice VI
    for iT = 1:nTbs
        if ~isempty(dindsVI(iT)) && dindsVI(iT)>0
            dzVI_m(iT) = R_Planet_m-R_sil_mean_m(dindsVI(iT)) - zVI_m(dindsVI(iT));
        elseif ~isempty(dindsV(iT)) && dindsV(iT)>0
            dzVI_m(iT) = 0;
            dzV_m(iT) = R_Planet_m-R_sil_mean_m(dindsV(iT)) - zVI_m(dindsV(iT));
        elseif ~isempty(dindsIII(iT)) && dindsIII(iT)>0
            dzV_m(iT) = 0;
            dzVI_m(iT) = 0;
            dzIII_m(iT) = R_Planet_m-R_sil_mean_m(dindsIII(iT)) - zVI_m(dindsIII(iT));
%         elseif ~isempty(dindsII(iT)) && dindsIII(iT)>0
%             dzV_m(iT) = 0;
%             dzVI_m(iT) = 0;
%             dzIII_m(iT) = R_Planet_m-R_sil_mean_m(dindsIII(iT)) - zVI_m(dindsIII(iT));
%             
        end
    end

    dzOcean_m(~dindsVI & ~dindsV) = R_Planet_m - R_sil_mean_m(~dindsVI & ~dindsV)-zI_m(~dindsVI & ~dindsV);

    zTotal_m = zI_m+dzOcean_m+dzIII_m+dzV_m+dzVI_m;

    disp(['Tb:                    ' num2str(Planet.Tb_K,'\t%0.0f')])
    disp(['z(km) ice I:           ' num2str(zI_m*1e-3,'\t%0.0f')])    
    disp(['z(km) ice III:         ' num2str(zIII_m*1e-3,'\t%0.0f')])    
    disp(['z(km) ice V:           ' num2str(zV_m*1e-3,'\t%0.0f')])    
    disp(['z(km) ice VI:          ' num2str(zVI_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) Ocean:          ' num2str(dzOcean_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) ice III:        ' num2str(dzIII_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) ice V:          ' num2str(dzV_m*1e-3,'\t%0.0f')])    
    disp(['dz(km) ice VI:         ' num2str(dzVI_m*1e-3,'\t%0.0f')])  
    disp(['dz(km) ice V + ice VI: ' num2str((dzV_m+dzVI_m)*1e-3,'\t%0.0f')])  
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
    [Q_Wm2(iT),deltaTBL_m(iT),eTBL_m(iT),Tc(iT),rhoIce(iT),alphaIce(iT),CpIce(iT),kIce(iT),CONVECTION_FLAG(iT)]=...
        ConvectionDeschampsSotin2001(Planet.Tsurf_K,Planet.Tb_K(iT),Pb_MPa(iT)/2,Zb2(iT),g_ms2(iT,1),2);
    if CONVECTION_FLAG(iT)
        %conductive upper layer
        inds_eTBL = find(z_m(iT,1:n_iceI)<=eTBL_m(iT));
        Pterm = P_MPa(iT,inds_eTBL)./P_MPa(iT,inds_eTBL(end));
        T_K(iT,inds_eTBL) = (Tc(iT).^(Pterm)).*(Planet.Tsurf_K.^(1-Pterm));
        %convective region
        for iconv = inds_eTBL(end)+1:n_iceI
            rho = 1000./getVspChoukroun2010(P_MPa(iT,iconv),T_K(iT,iconv-1),2);
            Cp = getCpIce(P_MPa(iT,iconv),T_K(iT,iconv-1),2);
            aK = 1.56e-4;
            T_K(iT,iconv) = T_K(iT,iconv-1)+aK*T_K(iT,iconv)/Cp/rho*deltaP*1e6;
        end
       % conductive lower layer
       inds_deltaTBL = find(z_m(iT,1:n_iceI)>=z_m(iT,n_iceI)-deltaTBL_m(iT));
        T_K(iT,inds_deltaTBL) = (Planet.Tb_K(iT).^(P_MPa(iT,inds_deltaTBL)./Pb_MPa(iT))).*(T_K(iT,inds_deltaTBL(1)-1).^(1-P_MPa(iT,inds_deltaTBL)./Pb_MPa(iT)));   
    end
end
% attenuation in ice
Hp_iceI = Seismic.g_aniso_iceI*T_K(iT,n_iceI);
QS_overfgamma_iceI = Seismic.B_aniso_iceI*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceI./T_K(:,1:n_iceI))/Seismic.LOW_ICE_Q;
Hp_iceIII = Seismic.g_aniso_iceIII*T_K(iT,n_iceI);
QS_overfgamma_iceIII = Seismic.B_aniso_iceIII*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceIII./T_K)/Seismic.LOW_ICE_Q;
Hp_iceV = Seismic.g_aniso_iceV*T_K;
QS_overfgamma_iceV = Seismic.B_aniso_iceV*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceV./T_K)/Seismic.LOW_ICE_Q;
Hp_iceVI = Seismic.g_aniso_iceVI*T_K(iT,n_iceI);
QS_overfgamma_iceVI = Seismic.B_aniso_iceVI*....
    exp(Seismic.gamma_aniso_iceI*Hp_iceVI./T_K)/Seismic.LOW_ICE_Q;

for iT = 1:nTbs % determine where the silicate should start
    indSil(iT) = find(R_sil_m(iT,:)==R_sil_mean_m(iT));
end
if nTbs>1 % this makes sure that the mantle length matrices are all the same length.   
    nsteps_mantle = [Params.nsteps_mantle Params.nsteps_mantle*ones(1,nTbs-1)+(indSil(1)-indSil(2:end))];% this needs to be reconsidered if nsteps!=100
else
    nsteps_mantle = Params.nsteps_mantle;
end

mprops = load(Seismic.mantleEOS);
rhofn = scatteredInterpolant([mprops(:,2),mprops(:,1)],mprops(:,3));
Cpfn = scatteredInterpolant([mprops(:,2),mprops(:,1)],mprops(:,4));
alphafn = scatteredInterpolant([mprops(:,2),mprops(:,1)],mprops(:,5));
VPfn = scatteredInterpolant([mprops(:,2),mprops(:,1)],mprops(:,6));
VSfn = scatteredInterpolant([mprops(:,2),mprops(:,1)],mprops(:,7));

%% Core
% Iron properties
    %gamma (fcc) Fe in Table 3 of C2006. derivatives are  rounded off values for alpha (bcc) Fe 
    % gamma Fe should be the correct phase between 1394°C (iron allotropes
    % wiki; 1667 K) and 912 °C (1185 K).  Minimum Tcore in C2006 is
    % 1200K.    
    rhoo_iron_kgm3 = 8e3;
    alpha_iron_K = 5;
%     Bulk modulus and derivatives?
    Kso_iron_GPa = 156;
    dKsdP_iron = 5;
    dKsdT_PaK = -0.040;
%     shear modulus and derivatives
    Go_iron_GPa = 76.5;
    dGdP_iron = 2;
    dGdT_iron_PaK = -0.023;

[Pcore_MPa,Tcore_K,rhocore_kgm3,M_above_core,VP_core_kms,VS_core_kms,g_ms2_core,Ks_iron_GPa,G_iron_GPa] = deal(zeros(nTbs,Params.nsteps_core));
for iT = 1:nTbs
    [Pmantle_MPa,rhomantle,M_above_mantle,VP_mantle_kms,VS_mantle_kms,g_ms2_sil] = deal(zeros(1,nsteps_mantle(iT)));
    Tmantle_K = T_K(iT,indSil(iT));
    r_mantle_m = linspace(R_sil_mean_m(iT),R_Fe_mean_m(iT),nsteps_mantle(iT));
    % cold case 
    if Planet.FeCore
        Tmantle_K = conductiveMantleTemperature(r_mantle_m,R_Fe_mean_m(iT),R_sil_mean_m(iT),Planet.kr_mantle,Planet.rho_sil_withcore_kgm3,Tmantle_K,Planet.Qmantle+Planet.QHmantle,0);
    else
        Tmantle_K = conductiveMantleTemperature(r_mantle_m,0,R_sil_mean_m(iT),Planet.kr_mantle,rho_sil_kgm3(iT,C2mean(iT)),Tmantle_K,Planet.Qmantle+Planet.QHmantle,0);
    end
    % COMMON FUNCTIONALITY HERE IS TO PROPAGATE PRESSURE, GRAVITY AND
    % DENSITY DOWNWARD FOR A CONDUCTIVELY COOLING LAYER WITH NO SOLID STATE
    % CONVECTION
    g_ms2_sil(1) = g_ms2(iT,C2mean(iT));
    M_above_mantle(1) = M_above_kg(iT,C2mean(iT));
    Pmantle_MPa(1) = P_MPa(iT,C2mean(iT)); 
    rhomantle(1) = 1e-3*rhofn(Pmantle_MPa(1)*10,Tmantle_K(1));%convert to g/mL
    for ij = 2:nsteps_mantle(iT)
         g_ms2_sil(ij)=(Gg*(M_Planet_kg-M_above_mantle(ij-1))/r_mantle_m(ij-1)^2);
         Pmantle_MPa(ij) = Pmantle_MPa(ij-1)+1e-6*(1000*rhomantle(ij-1))*g_ms2_sil(ij)*(r_mantle_m(ij-1)-r_mantle_m(ij));
         M_above_mantle(ij) = M_above_mantle(ij-1)+4/3*pi*(r_mantle_m(ij-1)^3-r_mantle_m(ij)^3)*1e3*rhomantle(ij-1);
         rhomantle(ij) = 1e-3*rhofn(Pmantle_MPa(ij-1)*10,Tmantle_K(ij-1)); %1e-3 convert to g/mL        
    end
    VP_mantle_kms = VPfn(Pmantle_MPa*10,Tmantle_K);
    VS_mantle_kms = VSfn(Pmantle_MPa*10,Tmantle_K);
    
%%  insert a low velocity layer in upper part of mantle
%     indsLow = find(r_mantle_m(1)-r_mantle_m<=5000);strLow = 'LowVUpperMantle5km';
%     indsLow = find(r_mantle_m(1)-r_mantle_m<=20000); strLow = 'LowVUpperMantle20km';
    indsLow = find(r_mantle_m(1)-r_mantle_m<=30000); strLow = 'LowVUpper30kmMantle';
    VP_mantle_kms(indsLow) = 4/7* VP_mantle_kms(indsLow);
    VS_mantle_kms(indsLow) = 4/7* VS_mantle_kms(indsLow);
    
    % AS ABOVE, NOW FOR THE METALLIC CORE       
    g_ms2_core(iT,1) = g_ms2_sil(nsteps_mantle(iT));
    M_above_core(iT,1) = M_above_mantle(nsteps_mantle(iT));
    r_core_m(iT,:) = linspace(R_Fe_mean_m(iT),0,Params.nsteps_core);
    Tcore_K(iT,:) = linspace(Tmantle_K(end),1.01*Tmantle_K(end),Params.nsteps_core);
    Pcore_MPa(iT,1) = Pmantle_MPa(nsteps_mantle(iT)); 
    Ks_iron_GPa(iT,1) = Kso_iron_GPa+1e-3*Pcore_MPa(iT,1)*dKsdP_iron+1e-9*Tcore_K(iT,1)*dKsdT_PaK;
    G_iron_GPa(iT,1) = Go_iron_GPa+1e-3*Pcore_MPa(iT,1)*dGdP_iron+1e-9*Tcore_K(iT,1)*dGdT_iron_PaK;
    rhocore_kgm3(iT,1) = rhoo_iron_kgm3; 
    for ij = 2:Params.nsteps_core
         g_ms2_core(iT,ij)=(Gg*(M_Planet_kg-M_above_core(iT,ij-1))/r_core_m(iT,ij-1)^2);
         Pcore_MPa(iT,ij) = Pcore_MPa(iT,ij-1)+1e-6*rhocore_kgm3(iT,ij-1)*g_ms2_core(iT,ij)*(r_core_m(iT,ij-1)-r_core_m(iT,ij));
         Ks_iron_GPa(iT,ij) = Kso_iron_GPa+1e-3*Pcore_MPa(iT,ij-1)*dKsdP_iron+1e-9*Tcore_K(iT,ij-1)*dKsdT_PaK;
         G_iron_GPa(iT,ij) = Go_iron_GPa+1e-3*Pcore_MPa(iT,ij-1)*dGdP_iron+1e-9*Tcore_K(iT,ij-1)*dGdT_iron_PaK;
         rhocore_kgm3(iT,ij) = rhocore_kgm3(iT,ij-1)*(1+1./(1e-3*Pcore_MPa(iT,ij)*Ks_iron_GPa(iT,ij))); 
         M_above_core(iT,ij) = M_above_core(iT,ij-1)+4/3*pi*(r_core_m(iT,ij-1)^3-r_core_m(iT,ij)^3)*rhocore_kgm3(iT,ij-1);
    end
    VS_core_kms(iT,:) = 1e-3*sqrt(G_iron_GPa(iT,:)*1e9./rhocore_kgm3(iT,:));
    VP_core_kms(iT,:) = 1e-3*sqrt(Ks_iron_GPa(iT,:)*1e9./rhocore_kgm3(iT,:)+4/3*(VS_core_kms(iT,:)*1e3).^2);    
    
    interior(iT).r_mantle_m = r_mantle_m;
    interior(iT).Dmantle = rhomantle;
    interior(iT).Tmantle_K = Tmantle_K;
    interior(iT).Pmantle_MPa = Pmantle_MPa;
    interior(iT).VS_mantle_kms = VS_mantle_kms;
    interior(iT).VP_mantle_kms = VP_mantle_kms;
        % compute seismic attenuation, as per C2006, 
    % Rb = 8.314462175; % J/K/mol
    Tm_K = 273.15+Tm_p_Hirschmann2000(1e-3*Pmantle_MPa); %input in GPa
    Hp_mantle = Seismic.g_aniso_mantle*Tm_K;
    interior(iT).QS_overfgamma = Seismic.B_aniso_mantle*exp(Seismic.gamma_aniso_mantle*Hp_mantle./Tmantle_K);
end

%% create the table of depths and heat flux values
% multirow{4}{*}{5 Wt\%} &Ice Ih &-&141& 127  & 96 & 63  & 24 \\
% & Liquid &-&63&92&128&244&539\\
% & Ice III &-&29&-&-&-&-\\
% & Ice V &-&158&134&49&-&-\\
% & Ice VI &-&409&411&411&355&237\\

Tb_str = getTableStr(Planet.Tb_K,Planet.Tb_K); %multiplying by 1e6 as a kluge to use the getDStr function
Qb_str = getTableStr(Planet.Tb_K,Qb*1e3); %multiplying by 1e6 as a kluge to use the getDStr function
Qc_str = getTableStr(Planet.Tb_K,Q_Wm2*1e3); %multiplying by 1e6 as a kluge to use the getDStr function
dzI_str = getTableStr(Planet.Tb_K,zI_m*1e-3);
dzOcean_str = getTableStr(Planet.Tb_K,dzOcean_m*1e-3);
dzIII_str = getTableStr(Planet.Tb_K,dzIII_m*1e-3);
dzV_str = getTableStr(Planet.Tb_K,dzV_m*1e-3);
dzVI_str = getTableStr(Planet.Tb_K,dzVI_m*1e-3);
R_Fe_str = getTableStr(Planet.Tb_K,R_Fe_mean_m*1e-3);
R_sil_str = getTableStr(Planet.Tb_K,(R_sil_mean_m)*1e-3);

disp('\hline')
disp(['\multirow{10}{*}{' num2str(wo) ' Wt\%}'])
if Planet.FeCore
    rhorockstr = ['\multicolumn{' num2str(nTbs) '}{c|}{' num2str(Planet.rho_sil_withcore_kgm3) '}'];
    disp(['&$\rho_{rock}$ (kg m$^{-3}$)&' rhorockstr '\\']);
else
    rhorockstr = getTableStr(Planet.Tb_K,rho_sil_kgm3(:,C2mean(:)));
    disp(['&$\rho_{rock}$ (kg m$^{-3}$)' rhorockstr{:} '\\']);
end
disp(['\cline{3-' num2str(nTbs+2) '}'])
disp(['&T$_{b}$ (K)  ' Tb_str{:} ' \\'])
disp(['&q$_{b}$ mW m$^{-2}$  ' Qb_str{:} ' \\'])
disp(['&q$_{c}$ mW m$^{-2}$  ' Qc_str{:} ' \\'])
disp(['&$D_{Ih}$ (km)' dzI_str{:} ' \\'])  
disp(['&$D_{ocean}$ (km)' dzOcean_str{:} ' \\'])
disp(['&$D_{III}$ (km)' dzIII_str{:} ' \\'])
disp(['&$D_{V}$ (km)' dzV_str{:} ' \\'])
disp(['&$D_{VI}$ (km)' dzVI_str{:} ' \\'])
disp(['&$R_{rock}$ (km)' R_sil_str{:} ' \\'])
if Planet.FeCore
    R_Fe_range_str = getTableStr(Planet.Tb_K,R_Fe_range_m*1e-3);
    R_sil_range_str = getTableStr(Planet.Tb_K,R_sil_range_m*1e-3);
    disp(['&R$_{core}$ (km)' R_Fe_str{:} ' \\'])
    disp(['&$\Delta$R$_{core} (km)$' R_Fe_range_str{:} ' \\'])
    disp(['&$\Delta$R$_{mantle}$ (km)' R_sil_range_str{:} ' \\'])
end
disp('\hline')



%% Sound Speeds in the Ice and Ocean
if Params.CALC_NEW_SOUNDSPEEDS
    disp('Computing sound speeds')
    velsIce = iceVelsGagnon1990(P_MPa,T_K);
    vfluid_kms = deal(zeros(nTbs,nsteps));
    for iT = 1:nTbs
        tic;
        ir = find(phase(iT,:)==0);
        vfluid_kms(iT,ir) = fluidSoundSpeeds(P_MPa(iT,ir),T_K(iT,ir),wo,Planet.Ocean.comp);
%         for ir = find(phase(iT,:)==0)
%             vfluid_kms(iT,ir) = fluidSoundSpeeds(P_MPa(iT,ir),T_K(iT,ir),wo,Planet.Ocean.comp);
%         end
        toc
        vfluid_kms(iT,vfluid_kms(iT,:)==0)=NaN;
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
    save([savefile 'Vels'],'velsIce','vfluid_kms');
else
    load([savefile 'Vels']);
end

%% Electrical Conductivity
if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
    LarionovMgSO4 = getLarionov1984MgSO4Conductivities;
    Pextrap = LarionovMgSO4.Pextrap_MPa;
    Textrap = LarionovMgSO4.Textrap_K;
    k_S_m_extrap = LarionovMgSO4.k_S_m_extrap_p01m;
    [Pgg,Tgg]=meshgrid(Pextrap,Textrap);
    k_S_mMgSO4p01Planet = NaN*ones(size(P_MPa)); 
    for iTb = 1:length(Planet.Tb_K)
        k_S_mMgSO4p01Planet(iTb,phase(iTb,:)==0) = interp2(Pextrap,Textrap,k_S_m_extrap,P_MPa(iTb,phase(iTb,:)==0),T_K(iTb,phase(iTb,:)==0),'spline');
    end
end

%% construct output for seismic modeling:
        % ice                           ocean                                            mantle         core 
[P_Planet_MPa,T_Planet_K,r_Planet_m,rho_pPlanet_kgm3,VP_Planet_kms,VS_Planet_kms,QS_overfgamma_Planet]=...
    deal(zeros(nTbs,indSil(1)-1+nsteps_mantle(1)+Params.nsteps_core));
%% Plot settings
LineWidth=1;
SoundLineWidth=1.5;
if Params.HOLD
    LineWidth = LineWidth+0.5;
    SoundLineWidth = SoundLineWidth+0.5;
    LineStyle='--';
else
    LineStyle='-';    
end
       
ymax = 1.05*R_Planet_m*1e-3;
nfig = 3000;
for iT = 1:nTbs
    % grab the indices for the ices, but only for the portion of the
    % interior calculation above the silicate interface
    H2Oinds = 1:indSil(iT)-1;
    indsLiquid = find(phase(iT,H2Oinds)==0);
    indsIII = find(phase(iT,H2Oinds)==3);
    indsV = find(phase(iT,H2Oinds)==5);
    indsVI = find(phase(iT,H2Oinds)==6);
    
    P_Planet_MPa(iT,:) = [P_MPa(iT,1:indSil(iT)-1) interior(iT).Pmantle_MPa Pcore_MPa(iT,:)];
    T_Planet_K(iT,:) = [T_K(iT,1:indSil(iT)-1) interior(iT).Tmantle_K Tcore_K(iT,:)];
    r_Planet_m(iT,:) = [r_m(iT,1:indSil(iT)-1)    interior(iT).r_mantle_m r_core_m(iT,:)];
    rho_pPlanet_kgm3(iT,:) = [rho_kgm3(iT,1:indSil(iT)-1) 1e3*interior(iT).Dmantle rhocore_kgm3(iT,:)];
    VP_Planet_kms(iT,:) = [velsIce.VIl_kms(iT,1:n_iceI) vfluid_kms(iT,indsLiquid) ...
        velsIce.VIIIl_kms(iT,indsIII) velsIce.VVl_kms(iT,indsV) velsIce.VVIl_kms(iT,indsVI) ...
        interior(iT).VP_mantle_kms VP_core_kms(iT,:)];
    VS_Planet_kms(iT,:) = [velsIce.VIt_kms(iT,1:n_iceI) 0*vfluid_kms(iT,indsLiquid) ...
        velsIce.VIIIt_kms(iT,indsIII) velsIce.VVt_kms(iT,indsV) velsIce.VVIt_kms(iT,indsVI) ...
        interior(iT).VS_mantle_kms VS_core_kms(iT,:)];
    QS_overfgamma_Planet(iT,:) = [QS_overfgamma_iceI(iT,:) 0*vfluid_kms(iT,indsLiquid)...
        QS_overfgamma_iceIII(iT,indsIII) QS_overfgamma_iceV(iT,indsV) QS_overfgamma_iceVI(iT,indsVI) ...
        interior(iT).QS_overfgamma Seismic.QScore*ones(1,Params.nsteps_core)];    
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        k_S_m(iT,:) = [0*QS_overfgamma_iceI(iT,:) k_S_mMgSO4p01Planet(iT,indsLiquid)...
            0*QS_overfgamma_iceIII(iT,indsIII) 0*QS_overfgamma_iceV(iT,indsV) 0*QS_overfgamma_iceVI(iT,indsVI) ...
            0*interior(iT).QS_overfgamma 0*Seismic.QScore*ones(1,Params.nsteps_core)];
        Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [P_Planet_MPa(iT,:)', T_Planet_K(iT,:)', r_Planet_m(iT,:)'*1e-3, rho_pPlanet_kgm3(iT,:)',VP_Planet_kms(iT,:)',VS_Planet_kms(iT,:)',QS_overfgamma_Planet(iT,:)' k_S_m(iT,:)'];
        header = sprintf('%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\t',...
            'P (MPa)','T (K)','r (km)','rho (kg m-3)','VP (km s-1)','VS (km s-1)','QS/gamma','k for 0.01 mol/kg MgSO4 (S m-1)');
    else
        Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma = [P_Planet_MPa(iT,:)', T_Planet_K(iT,:)', r_Planet_m(iT,:)'*1e-3, rho_pPlanet_kgm3(iT,:)',VP_Planet_kms(iT,:)',VS_Planet_kms(iT,:)',QS_overfgamma_Planet(iT,:)'];
        header = sprintf('%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t',...
            'P (MPa)','T (K)','r (km)','rho (kg m-3)','VP (km s-1)','VS (km s-1)','QS/gamma');
    end
    thissavestr = [savefile 'Zb' strLow num2str(1e-3*Zb2(iT),'%0.0f') 'km'];
%     save(thissavestr,'Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma','-ascii');
    dlmwrite([thissavestr '.txt'],header,'delimiter','');
     dlmwrite([thissavestr '.txt'],Wtpct_PMPaTKRkmRhokgm3VPkmsVSkmsQsoverfgamma,...
             'delimiter','\t',...
             'precision','%3.5e',...
             '-append');


%% plot the seismic data and attenuation
    figure(nfig+iT);
%     if ~Params.HOLD
%         clf
%     end
    set(gcf,'Position', [476   662   560   220],'Name',thissavestr)
    hp = subplot(1,3,1);
%     if Params.HOLD
%         hold on
%     end
    plot(VS_Planet_kms(iT,:)',r_Planet_m(iT,:)'*1e-3,...
        VP_Planet_kms(iT,:)',r_Planet_m(iT,:)'*1e-3,'--','LineWidth',LineWidth)
        set(gcf,'color','w')
        xlabel('Sound Speeds (km s^{-1})')
        ylabel(['r_{' Planet.name '} (km)']);
        set(gca,'xlim',[0 1.1*max(VP_Planet_kms(iT,:))],'ylim',[0 ymax])
        grid on

    
    hp = subplot(1,3,2);
%     if Params.HOLD
%         hold on
%     end
    plot(T_Planet_K(iT,:)',r_Planet_m(iT,:)'*1e-3,...
        rho_pPlanet_kgm3(iT,:)',r_Planet_m(iT,:)'*1e-3,'--','LineWidth',LineWidth);
    xlabel('Temperature (K), Density (kg m^{-3})');
    set(gca,'ylim',[0 ymax],'xlim',[0 max(rho_pPlanet_kgm3(iT,:))]);
    grid on
    
    subplot(1,3,3)
%     if Params.HOLD
%         hold on
%     end
    hp = plot(QS_overfgamma_Planet(iT,:)',r_Planet_m(iT,:)'*1e-3,'LineWidth',LineWidth);
    set(gca,'xscale','log','ylim',[0 ymax],'xlim',[10 1e7],'XTick',[10 100 1000 1e4 1e5 1e6 1e7])
    grid on
    xlabel('Q_S/\omega^{\gamma}')
    set(gcf,'color','w')
        
    saveas(gcf,['figures/' thissavestr 'QS'],Params.savefigformat);

end

%% create the legend that describes tb, z_b, and heat flux
lstr_2 = {};
for iT = 1:nTbs
    lstr_2 = {lstr_2{:} ['T_{b}:' num2str(Planet.Tb_K(iT),'%0.2f') ' K, q_{b}/q_{c}:'...
        num2str(1e3*Qb(iT),'%0.0f') '/' num2str(1e3*Q_Wm2(iT),'%0.0f') ' mW m^{-2}, z_{b}:' num2str(1e-3*Zb2(iT),'%0.0f') ' km']};
end

%%  plot profile with 4 subplots
%% 
if Params.foursubplots
Dsil_km=(R_Planet_m-R_sil_mean_m(:))*1e-3;

figure(229);
    ymaxScale = 1.01;
if ~Params.HOLD
    clf;
end
    set(gcf,'Position', [411    52   764   468])
    subplot(3,2,2);
%     subplot(1,2,2);
    hold on
for iT = 1:nTbs
    line(T_K(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineWidth',LineWidth,'LineStyle',LineStyle);
    hm = line(T_K(iT,indSil(iT)),z_m(iT,indSil(iT))*1e-3,'Color',Params.colororder(iT),'Marker','o');
    if ~Params.HOLD
        hm.MarkerFaceColor=Params.colororder(iT);
    end
end
if Params.HOLD
    ax = gca;
    if max(max(T_K))>ax.XLim
        ax.XLim(2)=max(max(T_K));        
    end
    if ymaxScale*max(Dsil_km)>ax.YLim
        ax.YLim(2) = ymaxScale*max(Dsil_km);
    end
else
    set(gca,'ydir','reverse','xlim',[245 max(max(T_K))],'ylim',[0 ymaxScale*max(Dsil_km)]);
end
xlabel('Temperature (K)');
ylabel('Depth (km)');
% for iT=1:nTbs
%     hline(Dsil_km(iT),Params.colororder(iT));
% end
box on;
    
    subplot(3,2,4); hold on;

for iT = 1:nTbs
    hp(iT) =  line(vfluid_kms(iT,Params.nsteps_iceI+1:indSil(iT)),z_m(iT,Params.nsteps_iceI+1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VIIIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVIl_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
    line(velsIce.VVIt_kms(iT,1:indSil(iT)),z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',SoundLineWidth);
end
if Params.HOLD
    ax = gca;
    if ymaxScale*max(Dsil_km)>ax.YLim
        ax.YLim(2) = ymaxScale*max(Dsil_km);
    end
else
    set(gca,'ydir','reverse','xlim',[1.3 5],'ylim',[0 ymaxScale*max(Dsil_km)]);%,'xlim',[1 4]);
end
ylabel('Depth (km)');
xlabel('Sound Speed (km s^{-1})');
box on



for iT=1:nTbs
    if ~isnan(velsIce.VVIl_kms(iT,indSil(iT)))
        velT = velsIce.VVIt_kms(iT,indSil(iT));
        velL = velsIce.VVIl_kms(iT,indSil(iT));
    elseif  ~isnan(velsIce.VVl_kms(iT,indSil(iT)))
        velT = velsIce.VVt_kms(iT,indSil(iT));
        velL = velsIce.VVl_kms(iT,indSil(iT));
    elseif  ~isnan(velsIce.VIIIl_kms(iT,indSil(iT)))
        velT = velsIce.VIIIt_kms(iT,indSil(iT));
        velL = velsIce.VIIIl_kms(iT,indSil(iT));
    elseif  ~isnan(velsIce.VIIl_kms(iT,indSil(iT)))
        velT = velsIce.VIIt_kms(iT,indSil(iT));
        velL = velsIce.VIIl_kms(iT,indSil(iT));
    elseif  ~isnan(vfluid_kms(iT,indSil(iT)))
        velT = vfluid_kms(iT,indSil(iT));
        velL = vfluid_kms(iT,indSil(iT));
    else
        velT = velsIce.VIt_kms(iT,indSil(iT));
        velL = velsIce.VIt_kms(iT,indSil(iT));
    end        
end

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
        subplot(3,2,6);hold on
        for iT = 1:nTbs
            line(k_S_mMgSO4p01Planet(iT,1:indSil(iT))*25,z_m(iT,1:indSil(iT))*1e-3,'Color',Params.colororder(iT),'LineStyle',LineStyle,'LineWidth',LineWidth);
        end
        if Params.HOLD
            ax = gca;
            if ymaxScale*max(Dsil_km)>ax.YLim
                ax.YLim(2) = ymaxScale*max(Dsil_km);
            end
        else
            set(gca,'ydir','reverse','xlim',[0.5 3],'ylim',[0 ymaxScale*max(Dsil_km)]);%,'xlim',[1 4]);
        end
        xlabel('Electrical Conductivity (S m^{-1} \times 25)');
        ylabel('Depth (km)');
        box on
    end



% subplot(1,2,1);
subplot(3,2,[1 3 5]);
hold all 
 if Params.HOLD
    ax = gca;
 end
R2ind = C2mean+npre; % map the index for R_sil_mean back to the indices for P, D, r_m, z_m, etc...
Psil_MPa = diag(P_MPa(:,R2ind(:)));
for iT=1:nTbs
    ht(iT)=  plot(P_MPa(iT,1:indSil(iT)),rho_kgm3(iT,1:indSil(iT)),Params.colororder(iT),'LineWidth',LineWidth,'LineStyle',LineStyle); 
    hm = plot(Psil_MPa(iT),interp1(P_MPa(iT,:),rho_kgm3(iT,:),Psil_MPa(iT)),[Params.colororder(iT) 'o']);
    if ~Params.HOLD
        hm.MarkerFaceColor=Params.colororder(iT);
    end
end

hw(1) = plot(Pref_MPa,rho_ref_kgm3(1,:),'k--');
hw(2) = plot(Pref_MPa,rho_ref_kgm3(2,:),'k--');
hw(3) = plot(Pref_MPa,rho_ref_kgm3(3,:),'k--');
hw(4) = plot(Pref_MPa,rho_ref_kgm3(4,:),'k--');


% hleg1 = legend([ht hw],{lstr_2{:},'0 wt%','5 wt%','10 wt%','15 wt%'});%,'20 wt%');
if Params.Legend
    hleg1 = legend([ht],lstr_2{:});
    set(hleg1,'location',Params.LegendPosition,'box','off')
end

box on;
 if Params.HOLD
    if ymaxScale*max(max((rho_kgm3(:,1:indSil(:)))))>ax.YLim
        ax.YLim(2) = ymaxScale*max(max((rho_kgm3(:,1:indSil(:)))));
    end
 else
    set(gca,'xlim',[0 ymaxScale*max(Psil_MPa)],'ylim',Params.ylim)
 end
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

    saveas(gcf,['figures/' savefile],Params.savefigformat);


%% plot profile with 2 subplots
else
figure(229);
ymaxScale = 1.01;
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
set(gca,'ydir','reverse','xlim',[0 5],'ylim',[0 ymaxScale*max(Dsil_km)]);%,'xlim',[1 4]);
    if Params.INCLUDE_ELECTRICAL_CONDUCTIVITY
        xlabel('Conductivity (S m^{-1} \times 25) and Sound Speed (km s^{-1})');
    else
        xlabel('Sound Speed (km s^{-1})');
    end
    ylabel('Depth (km)');
box on



for iT=1:nTbs
    if ~isnan(velsIce.VVIl_kms(iT,indSil(iT)))
        velT = velsIce.VVIt_kms(iT,indSil(iT));
        velL = velsIce.VVIl_kms(iT,indSil(iT));
    elseif  ~isnan(velsIce.VVl_kms(iT,indSil(iT)))
        velT = velsIce.VVt_kms(iT,indSil(iT));
        velL = velsIce.VVl_kms(iT,indSil(iT));
    elseif  ~isnan(velsIce.VIIIl_kms(iT,indSil(iT)))
        velT = velsIce.VIIIt_kms(iT,indSil(iT));
        velL = velsIce.VIIIl_kms(iT,indSil(iT));
    elseif  ~isnan(velsIce.VIIl_kms(iT,indSil(iT)))
        velT = velsIce.VIIt_kms(iT,indSil(iT));
        velL = velsIce.VIIl_kms(iT,indSil(iT));
    elseif  ~isnan(vfluid_kms(iT,indSil(iT)))
        velT = vfluid_kms(iT,indSil(iT));
        velL = vfluid_kms(iT,indSil(iT));
    else
        velT = velsIce.VIt_kms(iT,indSil(iT));
        velL = velsIce.VIt_kms(iT,indSil(iT));
    end        
end

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

set(gca,'ydir','reverse','xlim',[245 320],'ylim',[0 ymaxScale*max(Dsil_km)]);
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
set(gca,'xlim',[0 ymaxScale*max(Psil_MPa)])
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

    saveas(gcf,['figures/' savefile],Params.savefigformat);
end
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

%% properties of water ice
%Note: for VspChoukroun inds: liquid = 1, Ice I = 2, II = 3, III = 4, V = 5, VI = 6
function rho_kgm3 = getRhoIce(P_MPa,T_K,ind)
    % convert to appropriate phase using our nomenclature
    if ~(ind==5 || ind==6)
        ind = ind+1;
    end
    rho_kgm3 = 1000./getVspChoukroun2010(P_MPa,T_K,ind);
function Cp = getCpIce(P_MPa,T_K,ind)
    % convert to appropriate phase using our nomenclature
    if ~(ind==5 || ind==6)
        ind = ind+1;
    end
    Cp = CpH2O_Choukroun(P_MPa,T_K,ind);
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
        case 'NH3'
            phase = getIcePhaseNH3(P_MPa,T_K,w_pct);
    end
function Tfreeze_K = getTfreeze(P_MPa,wo,str_comp)
    switch str_comp
        case 'MgSO4'
%             if P_MPa<800
%                 W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
%                 mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from Wt%
%                 Tfreeze_K = TmeltMgSO4Bollengier(P_MPa,mo); % this is faster and more accurate
%             else
                Tfreeze_K = fzero(@(T) L_IceMgSO4(P_MPa,T,wo,1),[200 400]);
%             end
        case 'NH3'
            Tfreeze_K = fzero(@(T) L_IceNH3(P_MPa,T,wo,1),[200 350]);
    end
function Pfreeze_MPa = getPfreeze(T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'             
            Pfreeze_MPa = fzero(@(P) L_IceMgSO4(P,T_K,wo,1),[0 250]);
        case 'NH3'
            Pfreeze_MPa = fzero(@(P) L_IceNH3(P,T_K,wo,1),[0 250]);
    end

%% fluid properties
function [rho_kgm3,Cp,alpha_Km1]=fluidEOS(P_MPa,T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'
            W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
            mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from Wt%
            [rho_kgm3,Cp,alpha_Km1]=MgSO4_EOS2_planetary_smaller(mo,P_MPa,T_K-273.15);
        case 'NH3'
            [rho_kgm3,Cp,alpha_Km1,~]=NH3_EOS(P_MPa,T_K,wo);
    end
function vel_kms = fluidSoundSpeeds(P_MPa,T_K,wo,str_comp)
    switch str_comp
        case 'MgSO4'
            W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
            mo=1000./W_MgSO4./(1./(0.01.*wo)-1); %conversion to molality from Wt%
             vel_kms = MgSO4_EOS2_planetary_velocity_smaller_vector(mo,P_MPa,T_K-273.15); % assumes P and T are the same length.
%             vel_kms = MgSO4_EOS2_planetary_velocity_smaller(mo,P_MPa,T_K-273.15);
        case 'NH3'
            [~,~,~,vel_out]=NH3_EOS(P_MPa,T_K,wo);
    end

%% Core Size
function [C2MR2,R_fe] = CoreSize(rho_s,rho_fe,C_H2O,M_above_kg,R_s)
global M_Planet_kg R_Planet_m
try
    R_fe = fzero(@(R_fe) getR_fe(rho_s,rho_fe,C_H2O,M_above_kg,R_s,R_fe),[-1500e3 1500e3]);
catch
    R_fe = NaN;
end
C2MR2 = C_H2O+8/15*pi*((R_s.^5-R_fe.^5)*rho_s+rho_fe*R_fe.^5);%/M_Planet_kg/R_Planet_m^2;

    function zero_me = getR_fe(rho_s,rho_fe,C_H2O,M_above_kg,R_s,R_fe)
    global M_Planet_kg 
    drho = rho_fe-rho_s;
    zero_me = 4*pi/3*rho_fe*R_fe.^3 - (M_Planet_kg - M_above_kg-4*pi/3*rho_s*(R_s.^3-R_fe.^3));
%% formatting for output
function d_str = getTableStr(Tb,Xin)
d_str = {};
for iT = 1:length(Tb)
    if Xin(iT)
        d_str{iT} = [' &' num2str(Xin(iT),'%0.0f') ];
    else
        d_str{iT} = '& -';
    end
end