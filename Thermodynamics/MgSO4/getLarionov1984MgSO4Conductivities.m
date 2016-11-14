function LarionovMgSO4 = getLarionov1984MgSO4Conductivities


TwentyFive2ZeroDegC=0.525; % Hand and Chyba 2007
LarionovMgSO4.T_K = [273 298  323    348    373    398    423]; %addition of 273K as per Hand and Chyba 2007
LarionovMgSO4.P_MPa = [
      0.1
      98.1
      196.1
      294.2
      392.3
      490.3
      588.4
      686.5
      784.6];
  
mMgSO4p01.Mg = .01; mMgSO4p01.SO4 = .01;
mMgSO4p005.Mg = .005; mMgSO4p005.SO4 = .005;
str_EOS = 'Millero';
for iP = 1:9
    for iT = 1:5
        Larionov_p01m.rho_gmL(iP,iT) = get_aq_density(mMgSO4p01,0.1*LarionovMgSO4.P_MPa(iP),LarionovMgSO4.T_K(iT),str_EOS);
        Larionov_p005m.rho_gmL(iP,iT) = get_aq_density(mMgSO4p005,0.1*LarionovMgSO4.P_MPa(iP),LarionovMgSO4.T_K(iT),str_EOS);
    end
end
%kluge the densities by approximating that density at 373K is same at 398
%and 423.  This is a limitation of the density calculation.
 Larionov_p01m.rho_gmL(:,6) =  Larionov_p01m.rho_gmL(:,5);
 Larionov_p01m.rho_gmL(:,7) =  Larionov_p01m.rho_gmL(:,5);
 Larionov_p005m.rho_gmL(:,6) =  Larionov_p01m.rho_gmL(:,5);
 Larionov_p005m.rho_gmL(:,7) =  Larionov_p01m.rho_gmL(:,5);
        

Larionov_p01m.Equiv_cm2_Ohm_mol = [
                154.6   238.6   317.5   375.6   412.1   408.6
                169.7   252.5   334.2   404.0   452.3   472.5
                174.3   257.2   340.3   415.0   473.0   507.4
                173.8   256.3   340.8   417.1   480.7   525.0
                170.4   252.4   336.5   414.2   481.1   531.9
                164.5   245.9   329.4   407.2   475.8   530.8
                157.6   237.7   320.3   398.0   468.3   526.2   
                149.3   229.0   310.4   387.7   458.2   517.7
                139.8   219.5   229.9   376.3   446.4   508.4];
Larionov_p01m.k_Ohm_cm(:,1) = TwentyFive2ZeroDegC*Larionov_p01m.Equiv_cm2_Ohm_mol(:,1)*0.01/1000./Larionov_p01m.rho_gmL(:,1);
Larionov_p01m.k_Ohm_cm(:,2:7) = Larionov_p01m.Equiv_cm2_Ohm_mol*0.01/1000./Larionov_p01m.rho_gmL(:,2:end);
Larionov_p01m.k_S_m = 100*Larionov_p01m.k_Ohm_cm;


  Larionov_p005m.Equiv_cm2_Ohm_mol = [
                233.8   366.8   543.0   656.3   783.1   864.2
                241.7   368.9   506.8   646.1   777.0   860.5
                240.4   363.3   494.7   629.2   757.7   843.2
                234.0   354.0   479.1   609.1   733.4   823.0
                224.9   342.8   460.7   587.3   706.2   800.9
                214.0   330.5   442.4   564.0   679.5   779.2
                202.3   317.4   423.5   541.4   654.2   755.4
                190.6   304.1   404.6   518.5   628.1   730.1
                177.1   290.4   386.3   495.5   600.6   704.3];
Larionov_p005m.k_Ohm_cm(:,1) = TwentyFive2ZeroDegC*Larionov_p005m.Equiv_cm2_Ohm_mol(:,1)*0.005/1000./Larionov_p005m.rho_gmL(:,1);
Larionov_p005m.k_Ohm_cm(:,2:7) = Larionov_p005m.Equiv_cm2_Ohm_mol*0.005/1000./Larionov_p005m.rho_gmL(:,2:end);
Larionov_p005m.k_S_m = 100*Larionov_p005m.k_Ohm_cm;

LarionovMgSO4.Larionov_p005m = Larionov_p005m;
LarionovMgSO4.Larionov_p01m = Larionov_p01m;
% a helpful plot in the range of interest for icy satellites
colors = 'brkg';
% for iT = 1:3
% plot(LarionovMgSO4.P_MPa,LarionovMgSO4.Larionov_p01m.k_S_m(:,iT),[colors(iT) '-o'],LarionovMgSO4.P_MPa,LarionovMgSO4.Larionov_p005m.k_S_m(:,iT),[colors(iT) '--*'])
% end
figure(2296);clf
hold on
for iT = 1:3
plot(LarionovMgSO4.P_MPa,LarionovMgSO4.Larionov_p01m.k_S_m(:,iT),[colors(iT) '-o'],...
    LarionovMgSO4.P_MPa,LarionovMgSO4.Larionov_p005m.k_S_m(:,iT),[colors(iT) '--*']);
end
box on
xlabel('Pressure (MPa)')
ylabel('Electrical Conductivity (S m^{-1})')
legend('273 K; 0.01 g/kg','273 K; 0.005 g/kg','298 K; 0.01 g/kg','298 K; 0.005 g/kg','323 K; 0.01 g/kg','323 K; 0.005 g/kg')
