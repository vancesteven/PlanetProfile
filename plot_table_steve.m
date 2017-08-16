% script for plotting werami table
% dependencies:
%          function read_perplex_weramitable.m
%          function smooth2a.m
% note: the script is not general, but appositely build to reproduce figures of specific werami table
%       as shown in Fig. ... of ...
% 
% IMPORTANT: two things need to be changed internally: 1. table name (row 12-15), 2.location and model name (for P-T profile) (row 97-98)
clear all

% 1. read table epyro_ water saturated
%name='epyrohp_sat_1.tab';
%name='echonhp_sat_1.tab';
name = 'echonhp_sat_6GPa.tab';
% name='echon_1.tab';
%name='epyro_1.tab';

printo=input(' give 1 if you wanto to save pdf plots  ');
poro=input(' give 1 if you wanto to apply porosity corrections  ');
%ane=input(' give 1 if you wanto to apply anelasticity correction  '); % you can add here your function to apply Q correction: I can prepare for you if you want 
smoo=input(' give 1 if you wanto to apply smoothing  ');
prot=input(' give 1 if you wanto to show properties along a P-T profile (file lu_pt) ');

[header,out] = read_perplex_weramitable(name);
if (poro==1) % apply porosity corrections for density, VP and VS
    cout=porosity_correction(out);
    clear out;
    out=cout;
    clear cout;
end
out.t=out.t-273; % degrees in Celsius

% density P-T table
if smoo==1
den=smooth2a(out.den,5); % smoothing tables
vp=smooth2a(out.vp,5); % smoothing tables
vs=smooth2a(out.vs,5); % smoothing tables
% smooth table, besides looking better graphically, are also more
% consistent with the possible resolution that is possible to achieve from
% geophysical surface observations
else
den=out.den; % smoothing tables
vp=out.vp; % smoothing tables
vs=out.vs; % smoothing tables 
end

figure(1), hold on
otitle='Pyrolite - water saturated ';
title (['\fontsize{15}', otitle])
pcolor(out.p,out.t,den);
shading interp
kk=[1500 3500];
caxis(kk)
xlabel('\fontsize{14} Pressure (GPa)')
ylabel('\fontsize{14} Temperature (^oC)')
colormap(jet(20))
c=colorbar;
c.Label.String = '\fontsize{15} density (kg m^{-3})';
outnam=strcat('density_',name,'.pdf');
if (printo==1)
print('-dpdf','-r100','-painters',outnam)
end

figure(2), hold on
title (['\fontsize{15}', otitle])
pcolor(out.p,out.t,vp);
shading interp
kk=[2.5 8.2];
caxis(kk)
xlabel('\fontsize{14} Pressure (GPa)')
ylabel('\fontsize{14} Temperature (^oC)')
colormap(jet(57))
c=colorbar;
c.Label.String = '\fontsize{15} V_P (km s^{-1})';
outnam=strcat('vp_',name,'.pdf');
if (printo==1)
print('-dpdf','-r100','-painters',outnam);
end

figure(3), hold on
title (['\fontsize{15}', otitle])
pcolor(out.p,out.t,vs);
shading interp
kk=[0.3 4.6];
caxis(kk)
xlabel('\fontsize{14} Pressure (GPa)')
ylabel('\fontsize{14} Temperature (^oC)')
colormap(jet(43))
c=colorbar;
c.Label.String = '\fontsize{15} V_S (km s^{-1})';
outnam=strcat('vs_',name,'.pdf');
if (printo==1)
print('-dpdf','-r100','-painters',outnam)
end


if (prot==1)
    old=pwd;
%     cd /Users/fabio/Dropbox/ALL/WORK/PROCEDURES/PLANETS/FROM_STEVE_VANCE/MODELS/ % here give directory with planetary models
cd Europa
   [zp zt zr zrho zvp zvs zqs zkappa]=textread('EuropaProfile_MgSO40WtPctZb5km.txt','%f%f%f%f%f%f%f%f','headerlines',1); % load model
% in=load('lu_PT.dat'); % load lu_PT.dat
    pp=zp.*1e-3;tt=zt-273; % P in Gpa and degrees in Celsius
    cd (old)


pli=reshape(out.p,size(out.p,1)*size(out.p,2),1);
tli=reshape(out.t,size(out.t,1)*size(out.t,2),1);
dli=reshape(den,size(den,1)*size(den,2),1);
vpli=reshape(vp,size(vp,1)*size(vp,2),1);
vsli=reshape(vs,size(vs,1)*size(vs,2),1);

F=TriScatteredInterp(pli,tli,vpli);
% F=scatteredInterpolant(pli,tli,vpli);
vpin=F(pp,tt);
clear F
F=TriScatteredInterp(pli,tli,vsli);
vsin=F(pp,tt);
clear F
F=TriScatteredInterp(pli,tli,dli);
denin=F(pp,tt);
clear F
clear pli tli dli vpli vsli dli

pp=-pp; % dummy... to plot with increasing P toward bottom...
% properties: density, Vp VS 
figure(), hold on
title (['\fontsize{15}', otitle])
subplot(1,3,1), hold on
plot(denin,pp,'k');plot(denin,pp,'ok');
%plot(zrho,pp,'r');plot(zrho,pp,'or');

xlabel('\fontsize{15} density (kg m^{-3})');
ylabel('\fontsize{15} Pressure (GPa)');
box('on'); grid('on');
subplot(1,3,2), hold on
plot(vpin,pp,'k');plot(vpin,pp,'ok');
%plot(zvp,pp,'r');plot(zvp,pp,'or');
xlabel('\fontsize{15} V_P (km s^{-1})');
ylabel('\fontsize{15} Pressure (GPa)');
box('on'); grid('on');
subplot(1,3,3), hold on
plot(vsin,pp,'k');plot(vsin,pp,'ok');
%plot(zvs,pp,'r');plot(zvs,pp,'or');
xlabel('\fontsize{15} V_S (km s^{-1})');
ylabel('\fontsize{15} Pressure (GPa)');
box('on'); grid('on');

% T profile
figure(), hold on
title (['\fontsize{15}', otitle])
subplot(1,3,1), hold on
plot(tt,pp,'k');plot(tt,pp,'ok');
xlabel('\fontsize{15} Temperature (^o C');
ylabel('\fontsize{15} Pressure (GPa)');
box('on'); grid('on');

pp=-pp; %dummy... reverse again... sorry about this...
% T profile
figure(1);clf; hold on
plot(pp,tt,'k','LineWidth',3);
figure(2);clf; hold on
plot(pp,tt,'k','LineWidth',3);
figure(3);clf; hold on
plot(pp,tt,'k','LineWidth',3);

end