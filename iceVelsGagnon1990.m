function vels_kms = iceVelsGagnon1990(P_MPa,T_K,ind)

if nargin == 3
    [Vl_kms,Vt_kms] = getV(P_MPa,T_K,ind);
    vels_kms.Vl_kms = Vl_kms;
    vels_kms.Vt_kms = Vt_kms;
elseif nargin==2
    [VIl_kms,VIt_kms] = getV(P_MPa,T_K,2);
    [VIIl_kms,VIIt_kms] = getV(P_MPa,T_K,3);
    [VIIIl_kms,VIIIt_kms] = getV(P_MPa,T_K,4);    
    [VVl_kms,VVt_kms] = getV(P_MPa,T_K,5);
    [VVIl_kms,VVIt_kms] = getV(P_MPa,T_K,6);
    
    vels_kms = struct('VIl_kms',VIl_kms,'VIt_kms',VIt_kms,...
        'VIIl_kms',VIIl_kms,'VIIt_kms',VIIt_kms,...
        'VIIIl_kms',VIIIl_kms,'VIIIt_kms',VIIIt_kms,...
        'VVl_kms',VVl_kms,'VVt_kms',VVt_kms,...
        'VVIl_kms',VVIl_kms,'VVIt_kms',VVIt_kms);
else
    [VVl_kms,VVt_kms,VVIl_kms,VVIt_kms] = getvGagnon(P_MPa);
    vels_kms = struct('VVl_kms',VVl_kms,'VVt_kms',VVt_kms,...
        'VVIl_kms',VVIl_kms,'VVIt_kms',VVIt_kms);
end
    
function [Vl_kms,Vt_kms] = getV(P_MPa,T_K,ind)
%be careful with inds: liquid = 1 (not treated here), Ice I = 2, II = 3, III = 4, V = 5, VI = 6    
    rho_iceV = 1000./getVspChoukroun2010(P_MPa,T_K,ind);
    [Ks,mu]=getBulkShearMod(P_MPa,T_K,ind);
    Vt_ms = sqrt(mu*1e9./rho_iceV);    
    Vt_kms = Vt_ms*1e-3;
    Vl_kms = 1e-3*sqrt(Ks*1e9./rho_iceV+4/3*Vt_ms.^2);
        
%from Sotin1998 (compiling Shaw1986, Gagnon 1988, Polian1983) 
% GPa     oC      GPa     GPa
function [Ks,mu] = getBulkShearMod(P_MPa,T_K,ind)
%from Gagnon 1990, table II
Ksmu = [
    0      0
    9.96  3.6 % Ih
    13.89 6.2 % II
    9.87  4.6 % III
    14.19 6.1 % V
    18.14 7.5 % VI    
    ];
Ks = Ksmu(ind,1);
mu = Ksmu(ind,2);
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
function [VVl_kms,VVt_kms,VVIl_kms,VVIt_kms] = getvGagnon(P_MPa)
%from Gagnon et al. 1990
    VVl_kms = 4.108 +0.000203*P_MPa;
    VVt_kms = 2.198 +0.000006*P_MPa;
    VVIl_kms = 4.193 +0.000456*P_MPa;
    VVIt_kms = 1.987 +0.000463*P_MPa;
    