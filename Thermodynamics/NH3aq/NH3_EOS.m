function [rho_out,Cp_out,alpha_out,vel_out]=NH3_EOS(P_in,T_in,x_in)
 
%,Videal,Vex,Vw,VDH
%function [rho,Vs,vel,Cp,alpha,Videal,Vex]=NH3_EOS(P,T,x)
% usage:
%  [rho,Vs,Videal,Vex,Vw,vel,Cp,alpha]=NH3_EOS(P,T,x,parms)
%   units: P in MPa, T in oC, volumes in cc/g Cp in J/kg/C, rho in kg/m3
%   velocities in km/s
 
% set up the EOS on the parameter grid
%load AqAmmoniaEOS20140122T234600
load AqAmmoniaEOS20150612T211759 % 
warning('off','MATLAB:interp3:NaNstrip');

[Pc,Tc,xc]=meshgrid(P_in,T_in,x_in);
% [Pgrid,mgrid,Tgrid]=meshgrid(Pnodes,mnodes,Tnodes);
 
%VDH=interp3(Pgrid,mgrid,Tgrid,VDHg,Pc,mc,Tc,'spline');
%Vex=interp3(Pgrid,mgrid,Tgrid,Vexg,Pc,mc,Tc,'spline');
%Videal=interp3(Pgrid,mgrid,Tgrid,Videalg,Pc,mc,Tc,'spline');
%Vs=interp3(Pgrid,mgrid,Tgrid,Vsg,Pc,mc,Tc,'spline');
rho_out=interp3(Pgrid,mgrid,Tgrid,rhoc,Pc,Tc,xc,'spline');
vel_out=interp3(Pgrid,mgrid,Tgrid,velc,Pc,mc,Tc,'spline');
alpha_out=interp3(Pgrid,mgrid,Tgrid,alphac,Pc,Tc,xc,'spline');
Cp_out=interp3(Pgrid,mgrid,Tgrid,Cpc,Pc,Tc,xc,'spline');

