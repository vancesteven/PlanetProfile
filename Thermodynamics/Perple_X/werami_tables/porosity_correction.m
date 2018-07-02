function [cout] = porosity_correction(out,phi_surface)
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
   
% 2. get porosity as function of depth (Vitovtova et al. 2014) (modified) also implement Han et al. (put vit=2 for it)
if exist('phi_surface') 
    if phi_surface>0
        vit = 2;
    else
        vit = 1;
    end
else
    vit=1;
end
 if (vit==1) % (Vitovtova et al. 2014)
     disp('Applying Vitovtova et al. 2014 porosity model')
     por = NaN(size(dep,1),size(dep,2));
     limde=26; % below this depth decrease become linear until 80km
    for l = 1:size(dep,1)
        for m = 1:size(dep,2)
            if (dep(l,m)<=limde)
                por(l,m) = -0.65 -0.1*dep(l,m) + 0.0019*(dep(l,m)^2); % Vitovtova
                por(l,m) = 10^por(l,m);
            elseif ((dep(l,m)>limde)&&(dep(l,m)<=80))
                pora =   -0.65 -0.1*limde + 0.0019*(limde^2); % decrease linearly from 26 km
                pora = 10^pora;
                fac=(80-dep(l,m))/(80-limde); % 1 at 26 0 at 80 km
                por(l,m)=pora.*fac;
            elseif (dep(l,m)>80)
                por(l,m)=0;
            end
        end
    end
 elseif (vit==2)  % get porosity as function of depth (as in Han et al. 2014)
     disp(['Applying Han et al. 2014 porosity model with phi_surface = ' num2str(phi_surface)])
     por = NaN(size(out.p,1),size(out.p,2));
     con=6.15; % for sedimentary rocks...
     %con=1.5; % adjusted using Vitovnova et al. 
%      surface_por=input('give a value for surface porosity ');
    if nargin ==2
        surface_por = phi_surface;
    else
        surface_por=.8;
    end
     pc=0.25; % closure pressure (in GPa) this depends on material (can be higher according last experiments on crystalline rocks, saito et al. 2016)    
    for l = 1:size(dep,1)
        for m = 1:size(dep,2)
            por(l,m) = surface_por*exp(-(con*out.p(l,m))/pc);
        end
    end
 end
 
% 3. get permeability as function of depth (Kuang and Jiao 2014) (not used here, but maybe later useful for dynamical modeling)
     per = NaN(size(dep,1),size(dep,2),5);
%      From Table 1
%      logkr=-25.4; logks=-11.5; alph=0.25; % Crust in general
%      logkr=-25.4; logks=-11.5; alph=0.25; % Upper crust in general
%      logkr=-19.0; logks=-12.0; alph=1.8;  % Oceanic crust
%      logkr=-19.5; logks=-8.0; alph=0.45;  % Disturbed crust
%         for m = 1:size(dep,2)
%             per(l,m) = logkr + (logks-logkr) *(1+dep(l,m))^(-alph); % Kuang & Jiao
%             per(l,m) = 10^per(l,m);
%         end
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
     
     for l = 1:size(dep,1)
        for m = 1:size(dep,2)
            per(l,m,:) = logkr + (logks-logkr).*(1+dep(l,m)).^(-alph); % Kuang & Jiao
        end
     end
     per = power(10,per);
    cout.per = per;
% 4. compute average density using water filled pores (assuming total porosity equal to effective porosity)
          fluid_rho=1023; % kg/m^-3
     for l = 1:size(dep,1)
        for m = 1:size(dep,2)
          cout.den(l,m,:) = por(l,m)*fluid_rho +(1-por(l,m))*out.den(l,m); % average density  
        end
     end
     
% 5. correct VP and VS - use Wyllie (1958) for VP decrease and the double
% (in %) for VS (very approximate, partly based on Christensen 1989 data),
% they found an even stronger difference between VP and VS
% there are also recent data on various materials, Yu et al. 2016, but with void pore spaces.
           vpfl=1475; % in m/s
           vpi=out.vp.*1e3;% convert in m/s
           ivp = NaN(size(dep,1),size(dep,2));ivs = NaN(size(dep,1),size(dep,2));
       for l = 1:size(dep,1)
        for m = 1:size(dep,2)
          cout.vp(l,m)=(1./(por(l,m)./vpfl+(1-por(l,m))./vpi(l,m))).*1e-3; % in km/s
         end
       end
         clear vpi
        dvp=((cout.vp-out.vp)./out.vp).*100;  % percent variation of VP
        if min(dvp)>=-40
            incfa=2; % increase factor for VS, according to Christensen et al. 1989 (see also Christensen review 2004)
        else
            incfa = 1;
        end
        dvs=dvp.*incfa;
        cout.vs=out.vs.*(1+dvs./100);
        cout.por = por;
        
% use Ji et al. relations for Lame' paramters ?? work in progress
% it can be done, but model effects of P,T including porosity and rock properties... 
% I think that is better to keep them separate, but can be used to check
% our values...
%     vsi=out.vs.*1e3; vpi=out.vp.*1e3;% in m/s
%    mu=out.den.*vsi.^2; % (in Pa) from thermodynamic table
%    lamb=out.den.*vpi.^2-2.*mu; % from thermodynamic table
%     clear vsi vpi
%   [~,~,den,~,a,dldp,c,k] = textread('table_1_Ji10_sim.txt','%s%s%f%s%f%f%f%f','headerlines',1); % read table 1 Ji et al. 2010
%
%   mlam=a+(dldp).*out.p-c.*exp(-k.*out.p)-
%   
  
end

