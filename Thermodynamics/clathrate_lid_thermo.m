function [Q_converge, T_clath_ice , deltaTBL_m,eTBL_m,Tc]=clathrate_lid_thermo(Tsurf, Tbasal, P,clath_depth)

K_S_hc=[10000 15000 20000]; % clathrate lid
K_S_etbl=[15000 18500 20000]; % conductive lid thickness
%K_S_Tc=[256.5 258.4 259.5]-5;
K_S_Tc=[256.5 258.4 259.5]-2.3;
qs=[7.8 5.4 4.0]/10^3; % heat flux


Q_converge=interp1(K_S_hc,qs,clath_depth); % assume linear relationship from K&S 2020 to get heat flux
%eTBL_m=interp1(K_S_hc,K_S_etbl,clath_depth); %assume linear relation from K&S 2020 to get thermal boundary layer
Tc=interp1(K_S_hc,K_S_Tc,clath_depth);
kIce = getK_Andersson2005(P,Tc,'Ih','T');

T_clath_ice=Q_converge*clath_depth/0.5+Tsurf; % Fourier's law Q=-kDt/dz
ice_etbl=kIce*(Tc-T_clath_ice)./Q_converge; % assume Fourier's law to find thickness of ice shell to meet Tc given T_clath_ice interface; 
eTBL_m=clath_depth+ice_etbl;

% % %     
%Tc=T_clath_ice.*exp(Q_converge.*(eTBL_m-clath_depth)./632);%heat flux for a conductive lid Ojakangas and Stevenson, 1989 also appears in Barr2009)
% % %    
%Tc=interp1(K_S_hc,K_S_Tc,clath_depth);
 %kIce = getK_Andersson2005(P,Tc,'Ih','T');
 
 deltaTBL_m=kIce*(Tbasal-Tc)./Q_converge;
end

% input the surface temperature, temperature at ice ocean interface,
% Pressure, Ice shell thickness, clathrate depth and gravity)
% heat_flux clath for a clathrate lid. Loops over temperature range to find
% where the heat fluxes are allowed to converge
% % P1=P(n_clath)/2;
% % P2=P(n_clath+n_ice)/2;
% % 
% % for int=(floor(Tsurf-floor(Tsurf))+1:(Tbasal)-floor(Tsurf))
% %     T_clath_ice_check(int)=int+floor(Tsurf);
% %     
% %     Q_Wm2_clath(int)=0.5*(T_clath_ice_check(int)-Tsurf)/clath_depth;
% %    % [Q_Wm2_clath(int),deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I_clath(int)]=...
% %     %    ConvectionDeschampsSotin2001(Tsurf,T_clath_ice_check(int),P1,clath_depth,g,33);
% %     [Q_Wm2_ice(int),deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I_ice(int)]=...
% %         ConvectionDeschampsSotin2001(T_clath_ice_check(int),Tbasal,P2,h_ice-clath_depth,g,2);
% %     %dT=Q_Wm2(iT)/0.5*n_clath(iT);
% % end;
% % 
% % 
% % [Q_min min_ind]=min(abs(Q_Wm2_clath-Q_Wm2_ice));
% % 
% % % if convection_flag is zero, that means ice will not convect;
% % % issue is that there's jump in whether or not convection happens 
% % % should append; if clathrates are completely convecting 
% % 
% % if Q_min <(1/1000)% difference is within 0.1mW/m2
% %     T_clath_ice=T_clath_ice_check(min_ind);
% %     [Q_Wm2,deltaTBL_m,eTBL_m,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]=...
% %         ConvectionDeschampsSotin2001(T_clath_ice,Tbasal,P2,h_ice-clath_depth,g,2);
% %   if CONVECTION_FLAG_I_ice(min_ind)~=0
% %     eTBL_m=eTBL_m+clath_depth;
% %       Q_converge=Q_Wm2;
% %    CONVECTION_FLAG=1;
% %   else
% %     CONVECTION_FLAG=0;
% %     Tc = Tbasal-Tsurf;
% %     deltaTBL_m = 0;
% %     eTBL_m = 0;
% %     Q_converge=Q_Wm2_clath(min_ind);  
% %       
% %   end
% %     
% % else % represents a situation where the clathrates should be convecting and there is not a convergence in Q 
% %    % error('Current version is not designed to allow clathrates to mix in convection with ice. Reduce lid thickness or change to pure clathrates')
% %      
% %     
% %     count=0;
% %     % modifying Deschamps
% %     varstrs = {'water','Ih','II','III','V','VI'};
% %     Rg=8.314; % J/mol/K - Ideal Gas Constant
% %     E_Jmol = [0 60e3...
% %         mean([98 55])*1e3 mean([103 151])*1e3...
% %         136e3   110e3]; %assuming DS2001 value for Ice I, mean values for II and III, and high T value for Ice VI
% %      nu0 = [0 1e14 1e18 5e12 5e14 5e14]; % (viscosity at melting point) Pa s %ice I value independently selected by DS2001 is the one recommended in Fig. 4b of Durham et al. 1997.  Other values are also from that table, minimum, along the melting curve.  
% %     Dcond = [0 632 418 242 328 183]; % term for ice V was adapted from the scaling of D*T^-0.612 to D*T^-1 based on the Table 1 of Andersson and Inaba
% %     
% %     % c1 and c2 from numerical experiments of , Deschamps and Sotin (2000) as
% %     % summarized in DS2001, in connection with Eq. 9
% %     c1 = 1.43;
% %     c2 = -0.03;
% %     DeltaT = Tbasal-Tsurf; % assume full temperature profile
% %     Tlith = Tsurf + 0.3*(DeltaT); % approximation that scales 150K from DS2001 to the generic case.
% %     
% %     % try adding loop to find convergence;
% %        ci_ratio=clath_depth/h_ice; % percentage of ice shell that's clathrates. Should really be ratio in convective only
% %    % ci_ratio_orig=ci_ratio+1;
% %     %while abs(ci_ratio-ci_ratio_orig)>0.01;
% %        % ci_ratio_orig=ci_ratio;
% %        
% %         Bice = E_Jmol(2)/2/Rg/c1;
% %          E_Jmol(7)=90000; %Durham 2003
% %         Bclath = E_Jmol(7)/2/Rg/c1;
% %         
% %         B=Bice+ci_ratio*(Bclath-Bice);% get a B value between ice and clathrates based on relative thicknesses to scale Tc 
% %         C = c2*DeltaT;
% %         Tc = B*(sqrt(1+2/B*(Tbasal-C))-1); %DS2001 Eq. 18
% %         Aice = E_Jmol(2)/Rg/Tbasal; % dimensionless
% %         Aclath=E_Jmol(7)/Rg/Tbasal;
% %         A=Aice+ci_ratio*(Aclath-Aice);
% %         nuice = nu0(2)*exp(Aice*(Tbasal/Tc-1)); % DS2001 Eq. 11., also Durham et al. 1997 %viscosity
% %         nuclath=nu0(2)*20*exp(Aclath*(Tbasal/Tc-1));
% %          nu=nuice+ci_ratio*(nuclath-nuice);
% %         SF=SeaFreeze([P2,Tc],varstrs{2});
% %         
% %         rhoIce=SF.rho; % density
% %         alphaIce=SF.alpha; %coefficent of thermal expansion
% %         CpIce=SF.Cp; % specific heat
% %         kIce = getK_Andersson2005(P2,Tc,varstrs{2},'T'); % W/m/K
% %          SF=Helgerud_sI(P2,Tc);
% %         rhoclath=SF.rho;
% %         alphaclath=SF.alpha;
% %         Cpclath = polyval([3.19 2150],Tc); % values from Ning et al 2014
% %         kclath=0.5;% from Kalousova and Sotin 2020
% %          rho=rhoIce+ci_ratio*(rhoclath-rhoIce);
% %           alpha=alphaIce+ci_ratio*(alphaclath-alphaIce);
% %          k=kIce+ci_ratio*(kclath-kIce);
% %           Cp=CpIce+ci_ratio*(Cpclath-CpIce);
% %           
% % %          Q_crit=k*(Tc-Tsurf)/(h_ice);
% % %          TBL_m_crit=kIce*(Tbasal-Tc)/Q_crit; % only ice is used assuming its not a pure clathrate layer and there's at least TBL_m of ice
% % %         
% %          Kappa = k/rho/Cp; % W m2 s
% %     Ra=alpha*rho*g*DeltaT*(h_ice^3)/Kappa/nu; % DS2001 Eq. 4  %Kalosuova 2017 uses nu0. This uses a combo of ice and clathrate values
% %     Ra_crit=10^4;
% %     if Ra>Ra_crit % Convection % this may be a kluge; check nu and Kappa, and read the literature to confirm that 10^5 is considered sufficient as indicated by DS2001Fig3
% %         
% %         CONVECTION_FLAG_I=1;
% %         Ra_del = 0.28*Ra^0.21; % DS2001 Eq. 8
% %         KappaIce=kIce/rhoIce/CpIce; % W m2 s
% %         %deltaTBL_m=(nu*Kappa/alpha/rho/g/(Tbasal-Tc)*Ra_del)^(1/3); % thermal boundary layer thickness, DS2001 Eq. 19; represents bottom layer which is pure ice
% %        % Q_Wm2=k*(Tbasal-Tc)/deltaTBL_m; % DS2001 Eq. 20 % heat flux with the ocean using pure ice values 
% %         deltaTBL_m=(nuice*KappaIce/alphaIce/rhoIce/g/(Tbasal-Tc)*Ra_del)^(1/3); % thermal boundary layer thickness, DS2001 Eq. 19; represents bottom layer which is pure ice
% %         Q_Wm2=kIce*(Tbasal-Tc)/deltaTBL_m; % DS2001 Eq. 20 % heat flux with the ocean using pure ice values 
% %         %Q_Wm2=k*(Tbasal-Tc)/deltaTBL_m; % using mixture of k value
% %        T_clath_ice=(Q_Wm2*clath_depth/0.5)+Tsurf;
% %        %T_basal_conductive=Tbasal-Q_Wm2*deltaTBL_m/kIce;
% %        eTBL_m_ice=kIce*(Tc-T_clath_ice)/Q_Wm2;
% %        if eTBL_m_ice>0
% %        eTBL_m=clath_depth+ eTBL_m_ice;
% %        else
% %            error('Current version is not designed to allow clathrates to mix in convection with ice. Reduce lid thickness or change to pure clathrates')
% %      
% %        end
% %       % eTBL_m =0.5*(Tc-Tsurf)/Q_Wm2; % Eq. 21 % uses clathrate value for top ice shell with Fourier's law
% % %         if eTBL_m >= clath_depth % some ice is including in the conduction layer
% % %             ci_ratio_cond=clath_depth/(eTBL_m);
% % %              k=kIce+ci_ratio_cond*(kclath-kIce);
% % %             eTBL_m =k*(Tc-Tsurf)/Q_Wm2; 
% % %         else
% % %             error('Current version is not designed to allow clathrates to mix in convection with ice. Reduce lid thickness or change to pure clathrates')
% % %         end
% % %         
% % %         %deltaTBL_m=(nuice*KappaIce/alphaIce/rhoIce/g/(Tbasal-Tc)*Ra_del)^(1/3); % thermal boundary layer thickness, DS2001 Eq. 19; represents bottom layer which is pure ice
% % %         %Q_Wm2=kIce*(Tbasal-Tc)/deltaTBL_m; % DS2001 Eq. 20 % heat flux with the ocean using pure ice values 
% % %         %eTBL_m =0.5*(Tc-Tsurf)/Q_Wm2; % Eq. 21 % uses clathrate value for top ice shell with Fourier's law
% % %         if eTBL_m+deltaTBL_m>h_ice % for now verify the conductive lid thickness is less than thickness of ice shell
% % %             CONVECTION_FLAG_I=0;
% % %             Tc = DeltaT;
% % %             deltaTBL_m = 0;
% % %             eTBL_m = 0;
% % %             
% % %                 Q_Wm2=k*DeltaT/h_ice; % clathrates don't have temp dependent Kice so Fourier's law can be used.
% % %             
% % %         end
% % %         %   eTBL_m = kIce*(Tlith-Ttop)/Q_Wm2; % Eq. 22
% % %     else % Conduction
% % %         CONVECTION_FLAG_I=0;
% % %         Tc = DeltaT;
% % %         deltaTBL_m = 0;
% % %         eTBL_m = 0;
% % %         if ind<7
% % %             Q_Wm2=Dcond(ind)*log(Tbasal/Ttop)/h_ice; %heat flux for a conductive lid Ojakangas and Stevenson, 1989 also appears in Barr2009)
% % %         else
% % %             Q_Wm2=kIce*DeltaT/h_ice; % clathrates don't have temp dependent Kice so Fourier's law can be used.
% % %         end
% % %     end
% % %     if clath_depth<eTBL_m % if the clathrates are self contained within conductive
% % %      T_clath_ice=Q_Wm2*clath_depth/0.5+Tsurf; % T=Q*z/k+Tsurf if conductive
% % %     else
% % %         T_clath_ice=Tc; % assume convective temperature
% % %     end
% %     
% %     
% %     
% %     
% %     Q_converge=Q_Wm2;
% %     Ice_thickness=h_ice;
% % %      ci_ratio=(clath_depth-eTBL_m)/h_ice; %ratio of convecting clathrates to convecting ice
% % %      if ci_ratio<0; % all clathrates are in conductive layer
% % %          ci_ratio=0;
% % %      end
% %     end
% %     
% %     
% % %     T_clath_ice=T_clath_ice_check(min_ind);
% % %     [Q_Wm2_ice(int),deltaTBL_m,eTBL_m_test,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]=...
% % %                 ConvectionDeschampsSotin2001(T_clath_ice,Tbasal,P2,h_ice-clath_depth,g,2);
% % %     while Q_min > (1/1000) & count <100 %& eTBL_m_test <0 % if heat fluxes don't converge
% % %        count=count+1;
% % %         h_ice=h_ice+1000; % increase ice shell thickness by 500 meter increments
% % %         %Q_Wm2_clath=[];
% % %         Q_Wm2_ice=[];
% % %         for int=(floor(Tsurf-floor(Tsurf))+1:(Tbasal)-floor(Tsurf))
% % %             T_clath_ice_check(int)=int+floor(Tsurf)-1;
% % %             %Q_Wm2_clath(int)=0.5*(T_clath_ice_check(int)-Tsurf)/clath_depth;
% % %             % use T_clath_ice at "surf" temp, P2+count to accomadate extra
% % %             % pressure, h_ice-clath for thickness
% % %             [Q_Wm2_ice(int),deltaTBL_m,eTBL_m(int),Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG_I]=...
% % %                 ConvectionDeschampsSotin2001(T_clath_ice_check(int),Tbasal,P2+count,h_ice-clath_depth,g,2);
% % %             %dT=Q_Wm2(iT)/0.5*n_clath(iT);
% % %         end;
% % %         
% % %         
% % %         [Q_min min_ind]=min(abs(Q_Wm2_clath-Q_Wm2_ice));
% % %         eTBL_m_test=eTBL_m(min_ind);
% %    % end
% %     
% %     
% %     
% %     %
% % 
% % 
% % end