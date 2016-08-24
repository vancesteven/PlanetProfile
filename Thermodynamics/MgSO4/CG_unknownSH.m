function [G] = CG_unknownSH(i,S0,H0,P,T)

%%% Anchor points (melting points) from Choukroun and Grasset (2010)
%%% i(phase): 1(Ih), 2(II), 3(III), 4(liq), 5(V), 6(VI)
%%% dS0 is the change of entropy at the transition between the liquid and the ice
%%% T0 in K, P0 in MPa, dSO in J/mol/K, dH0 in J/mol
T0  = [251.16  252.32  256.16  273.15  256.16  273.31] ;
P0  = [209.90  300.00  350.10  0.1000  350.10  632.40] ;
dS0 = [-18.79  -21.82  -16.19  0.0000  -17.43  -18.78] ;

%%% Calculation of G 
dH = integral(@(T)CG_Cp(i,T),T0(i),T) + (integral(@(P)CG_V(i,P,T),P0(i),P).*10^3) ;
dS = integral(@(T)CG_CpT(i,T),T0(i),T) ;

H = H0 + dH ;
S = S0 + dS ;
G = H - T.*S ;
%   result from the V integral in kJ.kg-1 is mutliplied by 10^3 to match
%   the J.kg-1 of the other parts of the equation

% Conversion of the above result in kJ.kg-1 for plotting
H = H ./ 1000 ;
S = S ./ 1000 ;
G = G ./ 1000 ;

% Old test with the XSteam function now called in a separate script
% IAPWS = (XSteam('h_pT',(P.*10),(T-273.15))) - T .* (XSteam('s_pT',(P.*10),(T-273.15)))

end





function [V] = CG_V(i,P,T)

% Function to calculate the specific volumes of pure water and water ices
% following the approach of Choukroun and Grasset (2010)

% INPUT ARGUMENTS
% i is the considered phase: 1(Ih), 2(II), 3(III), 5(V), 6(VI) and 4(liq)
% P is the pressure in MPa
% T is the temperature in Kelvin

% OUTPUT ARGUMENT
% The function returns the specific volume in dm3.kg-1
% For reference, 1 m-3 equals 1 J.Pa-1 



% Parameters from Choukroun and Grasset (2010)
% Below are the actual values in the program, x for values that differs from the paper
% This set of parameters do not reproduce the Vs reported in the paper for ices II and VI
%            Ih       II x    III x    liq       V        VI x
%  V0   = [1.08600  0.85510  0.85460  0.81500  0.78300  0.74300] ;
%  Tref = [273.160  248.850  256.430  400.000  273.310  356.150] ;
%  a0   = [0.01900  0.02320  0.03750  0.10000  0.00500  0.02400] ;
%  a1   = [0.00750  0.10800  0.02030  0.00500  0.01000  0.02000] ;
%  b0   = [0.97400  0.99100  0.95100  1.00000  0.97700  0.96900] ;
%  b1   = [0.03020  0.07930  0.09700  0.28400  0.12000  0.05000] ;
%  b2   = [0.00395  0.00500  0.00102  0.00136  0.00160  0.00020] ;
% Below are the original values from the paper
 V0   = [1.08600  0.84250  0.85500  0.81500  0.78300  0.74300] ;
 Tref = [273.160  238.450  256.430  400.000  273.310  356.150] ;
 a0   = [0.01900  0.06000  0.03750  0.10000  0.00500  0.02400] ;
 a1   = [0.00750  0.00700  0.02030  0.00500  0.01000  0.00200] ;
 b0   = [0.97400  0.97600  0.95100  1.00000  0.97700  0.96900] ;
 b1   = [0.03020  0.04250  0.09700  0.28400  0.12000  0.05000] ;
 b2   = [0.00395  0.00220  0.00200  0.00136  0.00160  0.00102] ;


% Temperature dependance of the specific volume
zeta1 = 1 + a0(i)*tanh(a1(i)*(T-Tref(i))) ;

% Pressure dependance of the specific volume
zeta2 = b0(i) + b1(i)*(1-tanh(b2(i)*P)) ;

% Calculation of the specific volume
V = V0(i) .* zeta1 .* zeta2 ;

%old linear test for ice VI, works well actually.....
%V = 0.755 - P*0.00002 ;

end





function [Cp] = CG_Cp(i,T)

% Function to calculate the heat capacities of pure water and water ices
% following the approach of Choukroun and Grasset (2010)

% INPUT ARGUMENTS
% i is the considered phase: 1(Ih), 2(II), 3(III), 5(V), 6(VI) and 4(liq)
% T is the temperature in Kelvin

% OUTPUT ARGUMENT
% The function returns the heat capacity in J.kg-1.K-1



% Parameters from Choukroun and Grasset (2010)
c0 = [74.11  2200.0  820.0  0  700.0  940.0] ;
c1 = [7.560  0.0000  7.000  0  7.560  5.500] ;

% Calculation of the heat capacity of the liquid phase
if i == 4
    
    % High-temperature region
    if T >= 231
        Cp = 4190 + 9*exp(-0.11*(T-281.6)) ;

    % Low-temperature region
    else
        Cp = 2142511.11 - 35312.772*T + 53.606*(T.^2) + 2.025691667*(T.^3) ...
            - 0.012166*(T.^4) + 2.33191*(10^-6)*(T.^5) + 1.36462*(10^-7)*(T.^6) ...
            - 2.68836*(10^-10)*(T.^7) ;
    end

% Calculation of the heat capacities of the ice phases
else
    Cp = c0(i) + c1(i)*T ;
end

end





function [CpT] = CG_CpT(i,T)

% Function to calculate Cp/T for pure water and water ices
% following the approach of Choukroun and Grasset (2010)

% INPUT ARGUMENTS
% i is the considered phase: 1(Ih), 2(II), 3(III), 5(V), 6(VI) and 4(liq)
% T is the temperature in Kelvin

% OUTPUT ARGUMENT
% The function returns Cp/T in J.kg-1.K-2



% Parameters from Choukroun and Grasset (2010)
c0 = [74.11  2200.0  820.0  0  700.0  940.0] ;
c1 = [7.560  0.0000  7.000  0  7.560  5.500] ;

% Calculation of the heat capacity of the liquid phase
if i == 4
    
    % High-temperature region
    if T >= 231
        Cp = 4190 + 9*exp(-0.11*(T-281.6)) ;
        CpT = Cp./T ;

    % Low-temperature region
    else
        Cp = 2142511.11 - 35312.772*T + 53.606*(T.^2) + 2.025691667*(T.^3) ...
            - 0.012166*(T.^4) + 2.33191*(10^-6)*(T.^5) + 1.36462*(10^-7)*(T.^6) ...
            - 2.68836*(10^-10)*(T.^7) ;
        CpT = Cp./T ;
    end

% Calculation of the heat capacities of the ice phases
else
    Cp = c0(i) + c1(i)*T ;
    CpT = Cp./T ;
end

end
