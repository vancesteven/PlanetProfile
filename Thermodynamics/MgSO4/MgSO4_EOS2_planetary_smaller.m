function [rho,Cp,alpha]=MgSO4_EOS2_planetary_smaller(m,P,T)
 
%,Videal,Vex,Vw,VDH
%function [rho,Vs,vel,Cp,alpha,Videal,Vex]=MgSO4_EOS2(m,P,T)
% usage:
%  [rho,Vs,Videal,Vex,Vw,vel,Cp,alpha]=MgSO4_EOS(m,P,T,parms)
%   units: P in MPa, T in �C, volumes in cc/g Cp in J/kg/C, rho in gm/cc
%   velocities in km/s
 
% set up the EOS on the parameter grid

load MgSO4EOS2_planetary_smaller_20121116
Cpg=Cp;
rhog=rho;
alphag=alpha;
%% interpolate the parameters onto the grid of the input
[Pgg,mgg,Tgg]= meshgrid(P_smaller_MPa,m_smaller_molal,T_smaller_C);
[Pc,mc,Tc]=meshgrid(P,m,T);

%VDH=interp3(Pgg,mgg,Tgg,VDHg,Pc,mc,Tc,'spline');
%Vex=interp3(Pgg,mgg,Tgg,Vexg,Pc,mc,Tc,'spline');
%Videal=interp3(Pgg,mgg,Tgg,Videalg,Pc,mc,Tc,'spline');
%Vs=interp3(Pgg,mgg,Tgg,Vsg,Pc,mc,Tc,'spline');
rho=interp3(Pgg,mgg,Tgg,rhog,Pc,mc,Tc,'spline');
% vel=interp3(Pgg,mgg,Tgg,velg,Pc,mc,Tc,'spline');
alpha=interp3(Pgg,mgg,Tgg,alphag,Pc,mc,Tc,'spline');
Cp=interp3(Pgg,mgg,Tgg,Cpg,Pc,mc,Tc,'spline');

