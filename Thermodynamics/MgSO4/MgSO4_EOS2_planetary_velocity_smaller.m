function [vel]=MgSO4_EOS2_planetary_velocity_smaller(m,P,T)
% specialized from of MgSO4EOS2:
% function [rho,Vs,vel,Cp,alpha,Videal,Vex]=MgSO4_EOS2(m,P,T)
% usage:
%  [rho,Vs,Videal,Vex,Vw,vel,Cp,alpha]=MgSO4_EOS(m,P,T,parms)
%   units: P in MPa, T in oC, volumes in cc/g Cp in J/kg/C, rho in gm/cc
%   velocities in km/s
 
% set up the EOS on the parameter grid
 load MgSO4_EOS_parms_2012_26_17_LT
%load MgSO4EOS2_planetary_smaller_20121116

% % interpolate the parameters onto the grid of the input
%  Vo=interp2(Tg_C,Pg_MPa,Vog,T(:)',P(:),'spline');
%  p1=interp2(Tg_C,Pg_MPa,p1g,T(:)',P(:),'spline');
%  p2=interp2(Tg_C,Pg_MPa,p2g,T(:)',P(:),'spline');
%  Vw=interp2(Tg_C,Pg_MPa,Vwg,T(:)',P(:),'spline');
[Pgg,mgg,Tgg]= meshgrid(Pg_MPa,mg,Tg_C);
[Pc,mc,Tc]=meshgrid(P,m,T);
 
%VDH=interp3(Pgg,mgg,Tgg,VDHg,Pc,mc,Tc,'spline');
%Vex=interp3(Pgg,mgg,Tgg,Vexg,Pc,mc,Tc,'spline');
%Videal=interp3(Pgg,mgg,Tgg,Videalg,Pc,mc,Tc,'spline');
%Vs=interp3(Pgg,mgg,Tgg,Vsg,Pc,mc,Tc,'spline');
% rho=interp3(Pgg,mgg,Tgg,rhog,Pc,mc,Tc,'spline');
vel=interp3(Pgg,mgg,Tgg,velg,Pc,mc,Tc,'spline');
% alpha=interp3(Pgg,mgg,Tgg,alphag,Pc,mc,Tc,'spline');
% Cp=interp3(Pgg,mgg,Tgg,Cpg,Pc,mc,Tc,'spline');

