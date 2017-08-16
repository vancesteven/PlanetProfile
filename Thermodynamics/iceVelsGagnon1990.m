function vels_kms = iceVelsGagnon1990(P_MPa,T_K,ind)

if nargin == 3
    [vels_kms.Vl_kms,vels_kms.Vt_kms] = getV(P_MPa,T_K,ind);
elseif nargin==2
    [vels_kms.VIl_kms,vels_kms.VIt_kms,vels_kms.KsI_GPa,vels_kms.GsI_GPa] = getV(P_MPa,T_K,2);
    [vels_kms.VIIl_kms,vels_kms.VIIt_kms,vels_kms.KsII_GPa,vels_kms.GsII_GPa] = getV(P_MPa,T_K,3);
    [vels_kms.VIIIl_kms,vels_kms.VIIIt_kms,vels_kms.KsIII_GPa,vels_kms.GsIII_GPa] = getV(P_MPa,T_K,4);    
    [vels_kms.VVl_kms,vels_kms.VVt_kms,vels_kms.KsV_GPa,vels_kms.GsV_GPa] = getV(P_MPa,T_K,5);
    [vels_kms.VVIl_kms,vels_kms.VVIt_kms,vels_kms.KsVI_GPa,vels_kms.GsVI_GPa] = getV(P_MPa,T_K,6);
    
%     vels_kms = struct('VIl_kms',VIl_kms,'VIt_kms',VIt_kms,...
%         'VIIl_kms',VIIl_kms,'VIIt_kms',VIIt_kms,...
%         'VIIIl_kms',VIIIl_kms,'VIIIt_kms',VIIIt_kms,...
%         'VVl_kms',VVl_kms,'VVt_kms',VVt_kms,...
%         'VVIl_kms',VVIl_kms,'VVIt_kms',VVIt_kms);
else
    [vels_kms.VIl_kms,vels_kms.VIt_kms,vels_kms.VIIl_kms,vels_kms.VIIt_kms,vels_kms.VIIIl_kms,vels_kms.VIIIt_kms,vels_kms.VVl_kms,vels_kms.VVt_kms,vels_kms.VVIl_kms,vels_kms.VVIt_kms] = getvGagnon(P_MPa);
end
%%========================

function [Vl_kms,Vt_kms,Ks,Gs] = getV(P_MPa,T_K,ind)
%be careful with inds: liquid = 1 (not treated here), Ice I = 2, II = 3, III = 4, V = 5, VI = 6    
    if ind ==1
        rho_iceV = w_rho_ice(T_K-gsw_T0,10*(P_MPa-gsw_P0/1e5));
    else
        rho_iceV = 1000./getVspChoukroun2010(P_MPa,T_K,ind);
    end
    [Ks,Gs]=iceKsGs(P_MPa,T_K,ind);

    %playing with pressure dependent mu
%     deltamu=.001*P_MPa;  
    deltaGs= 0;
    
    Vt_ms = sqrt((Gs+deltaGs)*1e9./rho_iceV);  
    Vt_kms = Vt_ms*1e-3;
    Vl_kms = 1e-3*sqrt(Ks*1e9./rho_iceV+4/3*Vt_ms.^2);
    
%%========================
function [Ks,Gs]=iceKsGs(P_MPa,T_K,ind)
%mu is used below synonymously with Gs, as mu is the term used by Shaw 1986

%from Gagnon 1990, table II
% Ksmu = [
%     0      0
%     9.96  3.6 % Ih
%     13.89 6.2 % II
%     9.87  4.6 % III
%     14.19 6.1 % V
%     18.14 7.5 % VI    
%     ];
% from Shaw 1986
% Ksmu = [
%     0      0
%     9.5  3.3 % Ih
%     13.89 6.2 % II
%     8.5  4.8 % III
%     13.2 6.0 % V
%     16.0 8.5 % VI    
%     ];

% svance 2016 11/23
Ksmu = [
    0      0
    9.5  3.3 % Ih
    13.89 5.15 % II
    8.9  2.7 % III
    11.8 5.7 % V
    14.6 5 % VI    
    ];

Ks = Ksmu(ind,1);
Gs = Ksmu(ind,2);

% from Shaw 1986, Table II, mup is inferred. corrections in pressure are
% due to the fact that the intercept values are relative to a high
% pressure.  Adopted intercept values from Gagnon 1990 are larger than
% those of Shaw 1986
% Kspmup= [
%   0     0
%   5.3   0 % Ih
%   0     0 % II
%   5.7-.21   2.3-.21 %III
%   5.2-.34   1.1-.34
%   6.6-.62   1.1-.62
% ];
% svance 2016 11/23
Kspmup= [
  0     0 %liquid
  5.3   0 % Ih
  1.6     3.5 % II
  3.65   6.55 %III
  4.8   0.9 %V
  4.1   3 % VI
];
if ind==2 %from Gagnon 1988
%     Ks=0.1*(92.3+5.55*P_MPa*1e-2-0.26*(P_MPa*1e-2).^2);
    Ks=0.1*(92.4+3.3*P_MPa*1e-2-0.26*(P_MPa*1e-2).^2);
    Gs=0.04*(89.7+5.37*P_MPa*1e-2-0.25*(P_MPa*1e-2).^2);
else
    Ks=Ks+Kspmup(ind,1)*P_MPa*0.001;
    Gs=Gs+Kspmup(ind,2)*P_MPa*0.001;
end
%from Sotin1998 (compiling Shaw1986, Gagnon 1988, Polian1983) 
% GPa     oC      GPa     GPa
% PTKsMuI[
%     0   -25     9.5     3.3     
%     .21 -25     10.6    3.3
%     .14 -35.5   9.96    3.6
%     0   0
%     .2 -20.3
%     ];
% PTKsMuII[
%     0.283 -35.5     13.89   6.2];
% PTKsMuIII[
%     0.21 -25    8.5     4.8
%     0.34 -25    9.3     5.1
%     0.276 -27.2 9.87    4.6   
%     0.2 -22.5
%     0.34 -17
%     ];

% these are within 7% of each other
% PTKsMuV[
%     0.34 -25    13.2    6.0
%     0.62 -25    14.7    6.3
%     0.48 -35.5  14.19   6.1
% %     0.39 -13.6
% %     0.59 -1.6
%     ];
% PTKsMuVI[
%     0.62 -25    14.9    8.3
%     0.80 -25    16      8.5
%     0.777 -35.5 18.14   7.5
%     0.59 -3.2
%     0.79 12.6
%     1    25
%     2    25
%     ];

%%========================
function [VIl_kms,VIt_kms,VIIl_kms,VIIt_kms,VIIIl_kms,VIIIt_kms,VVl_kms,VVt_kms,VVIl_kms,VVIt_kms] = getvGagnon(P_MPa)
%from Gagnon et al. 1988 (I) and 1990
    VIl_kms = 3.914 +0.000452*P_MPa-0.0031*(P_MPa*1e-2).^2;
    VIt_kms = 1.997 +0.000250*P_MPa-0.0011*(P_MPa*1e-2).^2;
    VIIl_kms = 4.323 +0.000060*P_MPa;
    VIIt_kms = 2.085 +0.000669*P_MPa;
    VIIIl_kms = 3.461 +0.000853*P_MPa;
    VIIIt_kms = 1.507 +0.001705*P_MPa;    
    VVl_kms = 4.108 +0.000203*P_MPa;
    VVt_kms = 2.198 +0.000006*P_MPa;
    VVIl_kms = 4.193 +0.000456*P_MPa;
    VVIt_kms = 1.987 +0.000463*P_MPa;

    