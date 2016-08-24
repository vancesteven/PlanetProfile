function [Tmelt_K,phase] = TmeltMgSO4Bollengier(P_MPa,m_molal)
%%% Routine to invert S and H for 'CG_unknownSG' function
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
phases = [1 3 5 6] ;
% molaltoplot = [0 .1 .5 1 1.5 2] ;
% wm = Molality2WtPercent(m_molal,120.3);


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

%%% DIRECT CHEMICAL POTENTIAL
MH2O = 15.9994 + 2*1.00794 ;
MMGSO4 = 24.305 + 32.066 + 4*15.9994 ;          

%Calculation of two bracketing molalities with constant dnH2O
%Given molality m implies 1kg H2O + m mole salt
%New molalities calculated for 1kg+/-1mole of H2O and same m mole of salt
mlow = (m_molal) / ((((1000/MH2O)+0.1)*MH2O)/1000) ;
mhigh = (m_molal) / ((((1000/MH2O)-0.1)*MH2O)/1000) ;


%%% Exploration phase by phase
for ii = phases
%     disp(['Working on phase: ' num2str(ii)])    
    %%% Reference: getting melting temperature Tm at current P from Simon-Glatzel equations   
    [Tm] = SG(ii,P_MPa) ;
    SGout(ii,1) = P_MPa ;
    SGout(ii,2) = Tm ;

    %%% Exploration in temperature: resetting temperature convergence param. at current P
    %%% T is the melting temperature that will result from the intersection of G surfaces
    %%% chi2 is to estimate the deviation between the intersection and the SG equations
    %%% Tmin and Tmax (in K) are the exploration limits to be used for all melting curves
    Tmelt_K = 1e10 ;
    chi2 = 0 ;
    Tmin = 253.2 ;
    Tmax = 373.1 ;

    %%% Exploration in temperature
    while Tmelt_K > 1e9
        %%% Calculating dG at Tmin, (Tmin+Tmax)/2, and Tmax
        for ik = 1:3
            Ttest = Tmin + ((Tmax-Tmin)/2) * (ik-1) ;

            %G can now be calculated for the above molalities
            Gmlow = fnval(spMgSO4raw,{mlow,P_MPa,Ttest}) ;
            Gmhigh = fnval(spMgSO4raw,{mhigh,P_MPa,Ttest}) ;
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
            dmu(ik) = CG_unknownSH(ii,S0(ii),H0(ii),P_MPa,Ttest) - muwater ;
            %%% Solution found if convergence criterion is reached
            if abs(dmu(ik)) < crit
                Tmelt_K = Ttest ;
            end
        end

        %%% Loop security and temperature domain reduction if convergence crit. not reached
        if Tmelt_K > 1e9
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

%     disp(['Melting point found for phase ' num2str(ii) ': '...
%         num2str(P_MPa) ' MPa, ' num2str(Tmelt_K) ' K (dT = ' num2str(Tmelt_K-Tm) ' K)'])

    Tout(ii) = Tmelt_K;
    varm(ii) = abs((Tm-Tmelt_K)^2);
end
    phase = find(varm == min(varm([1 3 5 6])));
    Tmelt_K = Tout(phase);

