function [rho,Cp] = REFPROPsurfaceAmmoniaWater(x,P_MPa,T_oC)

D = refpropm('D','T',323.15,'P',1e2,'water','ammonia',[0.9 0.1])
%      Density of a 10% ammonia/water solution at 100 kPa and 323.15 K.