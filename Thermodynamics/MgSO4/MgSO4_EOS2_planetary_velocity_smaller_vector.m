function [vel]=MgSO4_EOS2_planetary_velocity_smaller_vector(m,P,T)
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
lT = length(T);
lP = length(P);
lm = length(m);
if lP~=lT
    error('length of P and T must be the same')
end
if lm==1
    mvec = m*ones(1,lT);
else
    mvec = m;
end
% [Pgg,mgg,Tgg]= meshgrid(Pg_MPa,mg,Tg_C);
% [Pc,mc,Tc]=meshgrid(P,m,T);
 
%VDH=interp3(Pgg,mgg,Tgg,VDHg,Pc,mc,Tc,'spline');
%Vex=interp3(Pgg,mgg,Tgg,Vexg,Pc,mc,Tc,'spline');
%Videal=interp3(Pgg,mgg,Tgg,Videalg,Pc,mc,Tc,'spline');
%Vs=interp3(Pgg,mgg,Tgg,Vsg,Pc,mc,Tc,'spline');
% rho=interp3(Pg_MPa,mg,Tg_C,rhog,P,mvec,T,'spline'); %gcc
vel=interp3(Pg_MPa,mg,Tg_C,velg,P,mvec,T,'spline'); %kms
% Ks = 1e-3*rho.*vel.^2; %GPa
% alpha=interp3(Pgg,mgg,Tgg,alphag,Pc,mc,Tc,'spline');
% Cp=interp3(Pgg,mgg,Tgg,Cpg,Pc,mc,Tc,'spline');

