function [cout] = porosity_correction_ice(out,por)

% made adjustments from original porosity_correction. Takes porosity, vp,
% vs, and density from inputs to determine changes to velocity curve.

% correct density, vp and vs of mineral physics table for porosity effects
% 1. Get porosity as function of pressure using porosity-depth relationship by Vitovtova et al. 2014
% and converting depth to pressure using PREM pressure profile.
% 2. Correct Vp/Vs using yu et al. 2016.
% 3. calculate density for fluid-filled porosity (effective porosity)
% inputs: p in GPa, t in K (not currently used), vp and vs in km/s,
% outputs: p,t unchanged, den in kg/m3, vp and vs in km/s
  cout.p=out.p; cout.t=out.t; % output same as input

% 1. convert pressure to depth of Min. phys. table according to simple PREM pressure profile
   dep=PREM_PRESSO(out.p); % in km
   
%waterout=SeaFreeze([out.p',out.t'],'water2');
   
  % disp(['Applying Han et al. 2014 porosity model with phi_surface = ' num2str(phi_surface)])
     %por = NaN(size(out.p,1),size(out.p,2));
    % 3. get permeability as function of depth (Kuang and Jiao 2014) (not used here, but maybe later useful for dynamical modeling)

     per = NaN(size(dep,1),size(dep,2),5);
      kskralpha = [
         -11.5  -25.4   0.25 % Crust in general
         -11.5  -25.4   0.25 % Upper crust in general
         -17.0  -21.2   1.50 % Low permeability upper crust
         -8.0   -19.5   0.45 % Disturbed crust    
         -12.0  -19.0   1.80 % Oceanic crust
         ];
     logkr = kskralpha(:,2);
     logks = kskralpha(:,1);
     alph = kskralpha(:,3);
      per = power(10,per);
    cout.per = per;
% 4. compute average density using water filled pores (assuming total porosity equal to effective porosity)
         % values are higher than used in previous version
%fluid_rho=1200; % kg/m^-3 % uses seafreeze to calculate density at given pressures and temp
%fluid_rho=1023; 
fluid_rho=0;
%vpfl=1475;%%
vpfl=0;
%vpfl=2500; % in m/s % assumes pores are filled with water
           vpi=out.vp.*1e3;% convert in m/s
           ivp = NaN(size(dep,1),size(dep,2));ivs = NaN(size(dep,1),size(dep,2));
      
     
     con=6.15; % for sedimentary rocks...
     %con=1.5; % adjusted using Vitovnova et al. ** No google search for
     %vitovnova
%      surface_por=input('give a value for surface porosity ');

    try
        surface_por = out.phi_surface;
    catch
        surface_por=.5; % larger assumption for Titan
    end
    pc=0.26; % % updated for Ice %closure pressure (in GPa) this depends on material (can be higher according last experiments on crystalline rocks, saito et al. 2016)    
    for l = 1:size(dep,1)
        for m = 1:size(dep,2)
            
           % por(l,m) = surface_por*exp(-(con*out.p(l,m))/pc);
           per(l,m,:) = logkr + (logks-logkr).*(1+dep(l,m)).^(-alph); % Kuang & Jiao
           
           % from Johnson et al 2017
           n(l,m)=(10e15).*exp(50e3./8.314.*(1./(out.t(l,m))-1./out.Tb));
           
           %viscosity with depth= viscoity in Gpa*
           %exp(Q*1000/R[1/T(z)-1/Tb) Q=50kJ/mol
           % porosity with depth =porosity at surface (-tpor*P/viscosity)
           por(l,m) = surface_por*exp(-(3600*24*365*50e3)*out.p(l,m)*10^9./n(l,m));
           
           %cout.den(l,m,:) = por(l,m)*fluid_rho +(1-por(l,m))*out.den(l,m);
            %cout.vp(l,m)=(1./(por(l,m)./vpfl+(1-por(l,m))./vpi(l,m))).*1e-3; % in km/s
           % cout.vp(l,m)=(1./(por(l,m)./vpfl+(1-por(l,m))./vpi(l,m))).*1e-3;
           cout.vp(l,m)=(1-por(l,m)).*vpi(l,m).*1e-3; % in km/s
        end
    end
    dvp=((cout.vp-out.vp)./out.vp).*100;  
     if min(dvp)>=-40
            incfa=2; % increase factor for VS, according to Christensen et al. 1989 (see also Christensen review 2004)
        else
            incfa = 1;
        end
        dvs=dvp.*incfa;
        cout.vs=out.vs.*(1+dvs./100);
       % cout.por = por;
end