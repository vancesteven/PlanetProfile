function sigma_Sm = getSigmaMgSO4_Pan(P_MPa,T_K)
% sigma_Sm = getSigmaMgSO4_Pan(wo,P_MPa,T_K)
% Steve Vance, 20191019
% returns values of electrical condcutivity for 10 wt% MgSO4 (aq) from .
% Regression data from lab experiments reproduce the orginal measurements
% with average misfit of 12%, as shown in the published supplement to the
% paper 
% "Electrical Conductivity of Aqueous Magnesium Sulfate at High Pressure 
% and Low Temperature with Application to Ganymede's Subsurface Ocean"
% Yilong Pan, Wenjun Yong, Richard A. Secco, Geophysical Research 

W_Mg = 24.3050; W_SO4 =  96.0636; W_MgSO4 = W_Mg + W_SO4;
EVAL = 0;
if ~exist('P_MPa')
    EVAL = 1;
    P_MPa = 0.1:10:1400.1;     
    T_K = 245:295; 
end
lP = length(P_MPa);lT = length(T_K);
wo = 10;
mo=1000./W_MgSO4./(1./(0.01.*wo)-1); % molality from wt%

if lP==lT
    [rho_kgm3,rho0_kgm3] = deal(ones(1,lT));
    for iT = 1:lT
        [rho_kgm3(iT),Cp,alpha_Km1]=MgSO4_EOS2_planetary_smaller(mo,P_MPa(iT),T_K(iT)-273.15);
        [rho0_kgm3(iT),Cp,alpha_Km1]=MgSO4_EOS2_planetary_smaller(0,P_MPa(iT),T_K(iT)-273.15);
    end
        c_M = mo*rho0_kgm3./rho_kgm3; % molarity
        logsig = getlogsig(P_MPa,T_K,c_M);    
else
    [c_M,logsig] = deal(ones(lP,lT));
    for iT = 1:lT
        [rho_kgm3,Cp,alpha_Km1]=MgSO4_EOS2_planetary_smaller(mo,P_MPa,T_K(iT)-273.15);
        [rho0_kgm3,Cp,alpha_Km1]=MgSO4_EOS2_planetary_smaller(0,P_MPa,T_K(iT)-273.15);
        c_M(:,iT) = mo*rho0_kgm3./rho_kgm3;
        logsig(:,iT) = getlogsig(P_MPa,T_K(iT),c_M(:,iT)');    
    end
end
sigma_Sm = 10.^logsig;

if EVAL
    set(0,'defaultfigurecolor',[1 1 1])
    set(0,'defaultAxesFontSize',16)

    figure(5);clf % reproduce Fig 5a from Pan et al. 2020
    cmap  = colormap(jet);
    cmap = cmap(20:end-20,:);
    colormap(cmap);
    hp = pcolor(T_K,P_MPa,sigma_Sm);shading interp
    
    hb = colorbar;
    caxis([0.15 3.4]);
    hb.Limits = [0.15 3.4];
    hb.Ticks = [0.15 3.4];
    
    hold on
    clevels = [0.17 0.62 1.1 1.5 2.0 2.4 2.9 3.3]; % from Pan et al. Fig. 5.
    [C,H]=contour(T_K,P_MPa,sigma_Sm,clevels);clabel(C)
    H.LineColor = 'k';

    set(gca,'ydir','reverse')
    ylabel('Pressure (MPa)')
    xlabel('Temperature (K)')
    ylim([0 1200]);
    xticks([250 260 270 280 290])
    box on
    
    figure(4);clf; % build a plot to compare with Fig. 4b of Pan et al. (2020). The regression data have upward concavity, whereas the measurements have downward concavity. 
    hold on; l_str = {};
    for iT = 1:round(lT/10):lT
        plot(P_MPa,sigma_Sm(:,iT))
        l_str = {l_str{:} [num2str(T_K(iT)) ' K']};
    end
    xlabel('Pressure (MPa)')
    ylabel('Conductivity (S/m)');
    hl = legend(l_str);
    hl.Box = 'on';
    box on
end
end % getSigmaMgSO4_Pan

function logsig = getlogsig(P_MPa,T_K,c_M)
G0 = getG0(P_MPa,T_K);
logsig = -3.1605 + 940.931./T_K + 0.8986.*log10(c_M)+2.*log10(G0) - log10(5*G0+2695.*c_M.^0.5);
end % getlogsig
function G0 = getG0(P_MPa,T_K)
rho = getrho(P_MPa,T_K);
G0 = 1918.37 - 100.51.*rho - 825071./T_K + 95550686./T_K.^2;
end %getG0
function rho_gmL = getrho(P_MPa,T_K)
rho_gmL = MgSO4_EOS2_planetary_smaller(0.9231,P_MPa,T_K-273.15)/1000;
end %getrho