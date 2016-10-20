function [Vex,VDH] = get_Vexs(m_molal,P_bar,T_K)
% T_K in Kelvin, P_bar in bars.
% m_molal should be a structure containing the individual ionic
% concentrations.  This has been vetted extensively for MgSO4.  It should
% be tested for other ions
% z = [2 1 2 -2 -1];
% I = 0.5*m*z'.^2;


% Following Marion et al (2005). 

% Calculate Vex =
% Av*I/b*ln(1+b*I^0.5)+2*R*T_K*sum(sum(mc*ma*(Bcav+sum(mc*zc)*Ccav)))
% Krumgalz (1995) has two additional terms in the last summation.  Marion
% (2005) truncates them

b = 1.2; % (kg/mol)^0.5, from Krumgalz (2000) page 1125
% 
%     Vex2MgCl = 
%     Vex2NaCl = get_Vex2Pitz(I,m(2),m(5),1,P_bar,T_K,'Na','Cl');
%     Vex2CaCl = get_Vex2Pitz(I,m(3),m(5),1,P_bar,T_K,'Ca','Cl');
%     Vex2MgSO4 = get_Vex2Pitz(I,m(1),m(4),2,P_bar,T_K,'Mg','SO4');
%     Vex2NaSO4 = get_Vex2Pitz(I,m(2),m(4),1,P_bar,T_K,'Na','SO4');
% Vex2 =  Vex2MgCl + Vex2NaCl + Vex2CaCl + Vex2MgSO4 + Vex2NaSO4;
Vex2 = get_Vex2Pitz(m_molal,P_bar,T_K);

lP = length(P_bar);  lT = length(T_K);
% Av from (Ananthaswamy and Atkinson 1984) (P_bar=1:1000 bars) (T_K=273:298 K)
for hi = 1:length(m_molal)
    I = get_ionic_strength(m_molal(hi));
    for ij = 1:lP
        for jk = 1:lT        
%             Av(ij,jk) = 3.73387-0.0289662*T_K(jk)+1.29461e-3*P_bar(ij)-5.62291e-6*T_K(jk).*P_bar(ij)+7.62143e-5*T_K(jk).^2+4.09944e-8*P_bar(ij).^2;
            Av(ij,jk) = getAv_BP79(P_bar(ij),T_K(jk)-273.15);
            Vex(hi,ij,jk) = Av(ij,jk).*I/b.*log(1+b*I.^0.5)+ squeeze(Vex2(hi,ij,jk));
            VDH(hi,ij,jk) = Av(ij,jk).*I/b.*log(1+b*I.^0.5);
        end
    end
end