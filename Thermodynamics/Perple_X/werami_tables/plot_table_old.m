% script for plotting werami table
% dependencies:
%          function read_perplex_weramitable.m
%          function smooth2a.m
% note: the script is not general, but appositely build to reproduce figures of specific werami table
%       as shown in Fig. ... of ...
% 
clear all

% 1. read table epyro_ water saturated
name='epyrohp_sat_1.tab';
%name='echonhp_sat_1.tab';
%name='echon_1.tab';
%name='epyro_1.tab';

printo=input(' give 1 if you wanto to save pdf plots  ');
poro=input(' give 1 if you wanto to apply porosity corrections  ');

[header,out] = read_perplex_weramitable(name);
if (poro==1) % apply porosity corrections for density, VP and VS
    cout=porosity_correction(out);
    clear out;
    out=cout;
    clear cout;
end
out.t=out.t-273; % degrees in Celsius

% density P-T table
den=smooth2a(out.den,10); % smoothing tables
vp=smooth2a(out.vp,10); % smoothing tables
vs=smooth2a(out.vs,10); % smoothing tables
% smooth table, besides looking better graphically, are also more
% consistent with the possible resolution that is possible to achieve from
% geophysical surface observations

figure(), hold on
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

figure(), hold on
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

figure(), hold on
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


