function [rho,alpha]=MgSO4_EOS2_planetary_velocity_smaller(m,P,T)
 
%,Videal,Vex,Vw,VDH
%function [rho,Vs,vel,Cp,alpha,Videal,Vex]=MgSO4_EOS2(m,P,T)
% usage:
%  [rho,Vs,Videal,Vex,Vw,vel,Cp,alpha]=MgSO4_EOS(m,P,T,parms)
%   units: P in MPa, T in �C, volumes in cc/g Cp in J/kg/C, rho in gm/cc
%   velocities in km/s
 
% set up the EOS on the parameter grid
% load MgSO4_EOS_parms6_AC
% load MgSO4_EOS_parms6_HD
% load MgSO4_EOS_parms6_AC
%  load MgSO4_EOS_parms8_IAPWS2 
% load MgSO4_EOS_parms
% load MgSO4_EOS_parms8_Cp
% load MgSO4_EOS_parms_2012_relaxed
% load MgSO4_EOS_parms_2012
% load MgSO4_EOS_parms_2012_R3
%  load MgSO4_EOS_parms_2012_27
%  load MgSO4_EOS_parms_2012_27_SV % this one has phi as phi/m
%  load MgSO4_EOS_parms_2012_27Mill
% load MgSO4_EOS_parms_2012_22fine
% load MgSO4_EOS_parms_2012_24fine
% load MgSO4_EOS_parms_2012_24fine46
% load MgSO4_EOS_parms_2012_25_20_LT
% load MgSO4_EOS_parms_2012_25_50_LT
% load MgSO4_EOS_parms_2012_25_17_LT
% load MgSO4_EOS_parms_2012_25_17_LT2
% load MgSO4_EOS_parms_2012_26_17_LT
load MgSO4EOS2_planetary_smaller_20121116

% % interpolate the parameters onto the grid of the input
%  Vo=interp2(Tg_C,Pg_MPa,Vog,T(:)',P(:),'spline');
%  p1=interp2(Tg_C,Pg_MPa,p1g,T(:)',P(:),'spline');
%  p2=interp2(Tg_C,Pg_MPa,p2g,T(:)',P(:),'spline');
%  Vw=interp2(Tg_C,Pg_MPa,Vwg,T(:)',P(:),'spline');
[Pgg,mgg,Tgg]= meshgrid(P_smaller_MPa,m_smaller_molal,T_smaller_C);
[Pc,mc,Tc]=meshgrid(P,m,T);
 
%VDH=interp3(Pgg,mgg,Tgg,VDHg,Pc,mc,Tc,'spline');
%Vex=interp3(Pgg,mgg,Tgg,Vexg,Pc,mc,Tc,'spline');
%Videal=interp3(Pgg,mgg,Tgg,Videalg,Pc,mc,Tc,'spline');
%Vs=interp3(Pgg,mgg,Tgg,Vsg,Pc,mc,Tc,'spline');
rho=interp3(Pgg,mgg,Tgg,rhog,Pc,mc,Tc,'spline');
%vel=interp3(Pgg,mgg,Tgg,velg,Pc,mc,Tc,'spline');
alpha=interp3(Pgg,mgg,Tgg,alphag,Pc,mc,Tc,'spline');
% Cp=interp3(Pgg,mgg,Tgg,Cpg,Pc,mc,Tc,'spline');

