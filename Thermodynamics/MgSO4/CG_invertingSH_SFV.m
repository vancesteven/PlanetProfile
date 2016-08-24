%%% Routine to invert S and H for 'CG_unknownSG' function

clearvars ;

%%% Model for water to be used
%load('The_spline_4_MgSO4.mat') ;
%load('IAPWS_Spline_250_1250_0_10000.mat') ;
%load('H2O_spline_0_5000_0_400_ver3.mat') ;
load('TryAgainSpline.mat') ;

%%% Pressure grid parameters, in MPa
% Pmin = [0.1  0  255  0  496  610 ] ;
% Pmax = [195  0  360  0  625  800] ;
Pmin = [0.1  0  200  0  350  610 ] ;
Pmax = [210  0  360  0  625  800] ;
nbP = 40 ;
colors = 'brgkmc';

%%% Convergence criterion for the intersection between G liq. and G ices
crit = 1e-3 ;

%%% Phase(s) for which S and H will be inverted
phasetoplot = [1 3 5 6] ;
molaltoplot = [0 .1 .5 1 1.5 2] ;
wtoplot = Molality2WtPercent(molaltoplot,120.3);

%%% Values to use with the Choukroun and Grasset representation of phase 4
%%% and actual CG volume values for ices
%S0(1) = -1.484895590000000e+03 ;
%S0(2) = -1.550670870000000e+03 ;
%S0(3) = -1.305093380000000e+03 ;
%S0(4) =  0 ;
%S0(5) = -1.370023840000000e+03 ;
%S0(6) = -1.310081570000000e+03 ;
%H0(1) = -1.775508950000000e+05 ;
%H0(2) = -1.135320270000000e+05 ;
%H0(3) = -0.104265700000000e+05 ;
%H0(4) =  0 ;
%H0(5) = -0.269506495900000e+05 ;
%H0(6) =  2.1250813500000000e+05 ;

%%% Starting values for the S0 and H0 of all phases (liquid: phase 4)
%%% Values to use with the Brown / IAPWS representation of phase 4 (pure water spline)
%%% and actual CG volume values for ices
%S0(1) = -1.435662151410000e+03 ;
%S0(2) = -1.550212380000000e+03 ;
%S0(3) = -1.090814460500000e+03 ;
%S0(4) =  0 ;
%S0(5) = -1.340912917200000e+03 ;
%S0(6) = -1.400305610000000e+03 ;
%H0(1) = -1.639072521000000e+05 ;
%H0(2) = -1.126204430000000e+05 ;
%H0(3) =  0.447400212700600e+05 ;
%H0(4) =  0 ;
%H0(5) = -0.192521063999988e+05 ;
%H0(6) =  1.878090121000000e+05 ;

%%% Starting values for the S0 and H0 of all phases (liquid: phase 4)
%%% Values to use with the Brown representation of phase 4 (MgSO4 spline)
%%% and actual CG volume values for ices
%S0(1) = -1.435662151410000e+03 ;
%S0(2) = -1.550212380000000e+03 ;
%S0(3) = -1.090814460500000e+03 ;
%S0(4) =  0 ;
%S0(5) = -1.340912917200000e+03 ;
%S0(6) = -1.400305610000000e+03 ;
%H0(1) = -1.6195072521000000e+05 ;
%H0(2) = -1.126204430000000e+05 ;
%H0(3) =  0.447400212700600e+05 ;
%H0(4) =  0 ;
%H0(5) = -0.192521063999988e+05 ;
%H0(6) =  1.897890121000000e+05 ;

%%% Starting values for the S0 and H0 of all phases (liquid: phase 4)
%%% Values to use with the Brown representation of phase 4 (MgSO4 spline number 2)
%%% and 2010 CG volume parameters for ice
S0(1) = -1.450662151410000e+03 ;
S0(2) = -1.550212380000000e+03 ;
S0(3) = -1.700814460500000e+03 ;
S0(4) =  0 ;
S0(5) = -1.340912917200000e+03 ;
S0(6) = -1.460305610000000e+03 ;
H0(1) = -1.678507252100000e+05 ;
H0(2) = -1.126204430000000e+05 ;
H0(3) = -1.115400212700600e+05 ;
H0(4) =  0 ;
H0(5) = -0.192521063999988e+05 ;
H0(6) =  1.712890121000000e+05 ;


%%% Exploration phase by phase
for ii = phasetoplot
    disp(['Working on phase: ' num2str(ii)])

    for im = 1:numel(molaltoplot)
    disp([' '])
    disp(['Working at molality: ' num2str(molaltoplot(im))])
    
    %%% Exploration in pressure: defining melting pressure P from 'ymin' to 'ymax' in 'ny' steps
    P = [Pmin(ii):((Pmax(ii)-Pmin(ii))/(nbP-1)):Pmax(ii)] ;
    for ij = 1:numel(P)

        %%% Reference: getting melting temperature Tm at current P from Simon-Glatzel equations   
        [Tm] = SG(ii,P(ij)) ;
        SGout(ii,ij,1) = P(ij) ;
        SGout(ii,ij,2) = Tm ;
        
        %%% Exploration in temperature: resetting temperature convergence param. at current P
        %%% T is the melting temperature that will result from the intersection of G surfaces
        %%% chi2 is to estimate the deviation between the intersection and the SG equations
        %%% Tmin and Tmax (in K) are the exploration limits to be used for all melting curves
        T = 1e10 ;
        chi2 = 0 ;
        Tmin = 253.2 ;
        Tmax = 373.1 ;

        %%% Exploration in temperature
        while T > 1e9

            %%% Calculating dG at Tmin, (Tmin+Tmax)/2, and Tmax
            for ik = 1:3
                Ttest = Tmin + ((Tmax-Tmin)/2) * (ik-1) ;
                %Gliq = CG(4,P(ij),Ttest) ;
                %%% SPLINES
                %Gliq = fnval(sp,{molaltoplot(im),P(ij)/1000,Ttest/1000})*1000 ;
                %Gliq = fnval(sp,{P(ij),Ttest}) ;
                %Gliq = fnval(sp2,{P(ij),Ttest}) ;
                %Gliq = fnval(spMgSO4raw,{molaltoplot(im),P(ij),Ttest}) ;
                %%% Correction from "kg solution" to "kg water"
                %Gliq = Gliq / (1 - (24.305+32.066+4*15.9994)*10^-3*molaltoplot(im)) ;
                %%% Test soustraction plutot que fraction
                %Gliq1 = fnval(spMgSO4raw,{molaltoplot(im),P(ij),Ttest}) ;
                %Gliq2 = fnval(spMgSO4raw,{0,P(ij),Ttest}) ;
                %dGliq = Gliq2 - Gliq1 ;
                %Gliq = Gliq2 - dGliq * (1 - (24.305+32.066+4*15.9994)*10^-3*molaltoplot(im)) ;
                %%% TEST IDEAL SOLUTION, R from J/mol/K to kJ/kg/K
                %R = 8.3144621 * 10^-3 * (1000/(2*1.00794+15.9994)) ;
                %Gliq = Gliq + R * Ttest * log(1) ;
                %%% TEST IDEAL MIXING, R from J/mol/K to kJ/kg/K
                %R = 8.3144621 * 10^-3 * (1000/(2*1.00794+15.9994)) ;
                %x = 1-1e10 ;
                %y = 1 - x ;
                %Gmix = R * Ttest * (x*log(x)+y*(log(y))) ;
                %Gliq = Gliq + Gmix ;

                
                %%% DIRECT CHEMICAL POTENTIAL
                                
                MH2O = 15.9994 + 2*1.00794 ;
                MMGSO4 = 24.305 + 32.066 + 4*15.9994 ;          
                
                %Calculation of two bracketing molalitites with constant dnH2O
                %Given molality m implies 1kg H2O + m mole salt
                %New molalities calculated for 1kg+/-1mole of H2O and same m mole of salt
                mlow = (molaltoplot(im)/1) / ((((1000/MH2O)+0.1)*MH2O)/1000) ;
                mhigh = (molaltoplot(im)/1) / ((((1000/MH2O)-0.1)*MH2O)/1000) ;
                %G can now be calculated for the above molalities
                Gmlow = fnval(spMgSO4raw,{mlow,P(ij),Ttest}) ;
                Gmhigh = fnval(spMgSO4raw,{mhigh,P(ij),Ttest}) ;
                %Output G are given for 1 kg of solution
                %Output G must be brought back to their respective nb of mole of water
                %G must be divided by actual number of H2O moles and scaled to original number
                nH2Olow = ((1000 / (1000 + mlow*MMGSO4))*1000) / MH2O ;
                Gmlow = Gmlow / nH2Olow * ((1000/MH2O)+0.1) ;
                nH2Ohigh = ((1000 / (1000 + mhigh*MMGSO4))*1000) / MH2O ;
                Gmhigh = Gmhigh / nH2Ohigh * ((1000/MH2O)-0.1) ;
                %mu water is the derivative of G
                muwater = (Gmlow-Gmhigh)/0.2 * (1000/MH2O) ;
                
                
                
                
                
                %%% FINAL CALCULATION
                %dmu(ik) = CG_unknownSH(ii,S0(ii),H0(ii),P(ij),Ttest) - Gliq ; 
                dmu(ik) = CG_unknownSH(ii,S0(ii),H0(ii),P(ij),Ttest) - muwater ;
                %%% Solution found if convergence criterion is reached
                if abs(dmu(ik)) < crit
                    T = Ttest ;
                end
            end

            %%% Loop security and temperature domain reduction if convergence crit. not reached
            if T > 1e9
                if ((dmu(1).*dmu(2)).*(dmu(2).*dmu(3))) > 0
                    [dmu(1), dmu(2), dmu(3)] ;
%                     error('No solution in the selected PT range!') ;
%                      disp('out of bounds')
                     break
                elseif (dmu(1)*dmu(2)) < 0
                    Tmax = Tmax - 0.5*(Tmax-Tmin) ;
                else
                    Tmin = Tmin + 0.5*(Tmax-Tmin) ;
                end
            end

        end
    
        disp(['Melting point found for phase ' num2str(ii) ': '...
            num2str(P(ij)) ' MPa, ' num2str(T) ' K (dT = ' num2str(T-Tm) ' K)'])
    
        CGout(ii,im,ij,1) = P(ij) ;
        CGout(ii,im,ij,2) = T ;
    
        chi2 = chi2 + abs((Tm-T)^2) ;
    end
    end
end

% % %         Pmin = [0.1  0  255  0  496  610 ] ;
% % %         Pmax = [180  0  360  0  625  800] ;
% % for ii = phasetoplot
% for ii = [5 6]
% for im = 1:numel(molaltoplot)
% %     P = [Pmin(ii):((Pmax(ii)-Pmin(ii))/(nbP-1)):Pmax(ii)] ;
%     Pm = squeeze(CGout(ii,im,:,1));
%     parfor ij = 1:length(Pm)
%         disp(['phase: ' num2str(ii) ...
%                 '; m: ' num2str(molaltoplot(im)) ...
%                 '; P: ' num2str(Pm(ij))])
%         Tv2014(ii,ij,im) = fzero(@(T) L_Ice(Pm(ij),T,wtoplot(im),1),[240 600]);
%     end
% end
% end
% Pm = squeeze(CGout(:,:,:,1));
% save('Tv2014MgSO420160612','Pm','Tv2014');
load Tv2014MgSO420160612

figure(522);clf;hold on
figure(523);clf;hold on
for ii = phasetoplot
    figure(522);
    plot(SGout(ii,:,1),SGout(ii,:,2),'k','LineWidth',1)
    for im = 1:numel(molaltoplot)
        Tmcalc = squeeze(CGout(ii,im,:,2));
        Tmcalc(Tmcalc>1e4) = nan;
        Pm = squeeze(CGout(ii,im,:,1));
        figure(522);
        plot(Pm,Tmcalc,colors(im))
        plot(Pm,Tv2014(ii,:,im),[colors(im) '--']);
    figure(523);
        plot(Pm,pdiff(Tv2014(ii,:,im)',Tmcalc),colors(im))
    end
end

    figure(522)
set(gca,'ylim',[200 max(max(Tmcalc))])

xlabel('Pressure (MPa)')
ylabel('Temperature (K)')
title('Melting Temperatures')

    figure(523)
xlabel('Pressure (MPa)')
yl = ylabel('$100\frac{T_{V2014}-T_{m}}{T_{m}}$ (K)');
yl.Interpreter = 'latex';
title('Deviations')
