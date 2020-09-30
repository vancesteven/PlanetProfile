function [Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG]=...
    ConvectionDeschampsSotin2001(Ttop,Tm,Pmid_MPa,h_m,g_ms2,ind)
% determine solid state convection for ices 
% based on Deschamps and Sotin 2001, Thermal Convection in Icy Satellites, J. Geophys. Res. 106(E3):5107-5121 (DS2001)

%inds 1-6 are [water, Ice Ih, Ice II, Ice III, Ice V, Ice VI]
varstrs = {'water','Ih','II','III','V','VI'};
Rg=8.314; % J/mol/K - Ideal Gas Constant
E_Jmol = [0 60e3...
    mean([98 55])*1e3 mean([103 151])*1e3...
    136e3   110e3]; %assuming DS2001 value for Ice I, mean values for II and III, and high T value for Ice VI
nu0 = [0 5e13 1e18 5e12 5e14 5e14]; % Pa s %ice I value independently selected by DS2001 is the one recommended in Fig. 4b of Durham et al. 1997.  Other values are also from that table, minimum, along the melting curve.
Dcond = [0 632 418 242 328 183]; % term for ice V was adapted from the scaling of D*T^-0.612 to D*T^-1 based on the Table 1 of Andersson and Inaba


% c1 and c2 from numerical experiments of , Deschamps and Sotin (2000) as
% summarized in DS2001, in connection with Eq. 9
c1 = 1.43; 
c2 = -0.03;
DeltaT = Tm-Ttop;
Tlith = Ttop + 0.3*(DeltaT); % approximation that scales 150K from DS2001 to the generic case.

B = E_Jmol(ind)/2/Rg/c1;
C = c2*DeltaT;
Tc = B*(sqrt(1+2/B*(Tm-C))-1); %DS2001 Eq. 18

A = E_Jmol(ind)/Rg/Tm; % dimensionless
nu = nu0(ind)*exp(A*(Tm/Tc-1)); % DS2001 Eq. 11., also Durham et al. 1997

rhoIce=1000./getVspChoukroun2010(Pmid_MPa,Tc,ind); % kg/m3
alphaIce = log(1000./rhoIce)-log(getVspChoukroun2010(Pmid_MPa,Tc-1,ind)); %1/K
CpIce = CpH2O_Choukroun(Pmid_MPa,Tc,ind); % J/kg/K
kIce = getK_Andersson2005(Pmid_MPa,Tc,varstrs{ind},'T'); % W/m/K
Kappa = kIce/rhoIce/CpIce; % W m2 s
Ra=alphaIce*rhoIce*g_ms2*DeltaT*h_m^3/Kappa/nu; % DS2001 Eq. 4
if Ra>10^5 % Convection % this may be a kluge; check nu and Kappa, and read the literature to confirm that 10^5 is considered sufficient as indicated by DS2001Fig3
    CONVECTION_FLAG=1;
    Ra_del = 0.28*Ra^0.21; % DS2001 Eq. 8
    deltaTBL_m=(nu*Kappa/alphaIce/rhoIce/g_ms2/(Tm-Tc)*Ra_del)^(1/3); % thermal boundary layer thickness, DS2001 Eq. 19
    Q_Wm2=kIce*(Tm-Tc)/deltaTBL_m; % DS2001 Eq. 20
    eTBL_m = kIce*(Tc-Ttop)/Q_Wm2; % Eq. 21 
%     eTBL_m = kIce*(Tlith-Ttop)/Q_Wm2; % Eq. 22
else % Conduction 
    CONVECTION_FLAG=0;
    Tc = DeltaT;
    deltaTBL_m = 0;
    eTBL_m = 0;
    Q_Wm2=Dcond(ind)*log(Tm/Ttop)/h_m; %
end
