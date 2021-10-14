function T_l_K = IronRichLiquidusRuckrieman(p_GPa,x_s)
% Calculation of liquidus based on Ruckrieman et al. 2018 (A1)
% p_GPa = pressure
% x_s = wt% of sulfur
% mole fraction of Iron Sulfide in the mixture
x_mol = 2.74189.*x_s; % inferred conversion 2.74189 = (atomic mass Fe + S) / (atomic mass S)
% I think the proper conversion would be
% (atomic mass Fe) / (atomic mass S) * (x_s) / (1-x_s)
% but the sources don't have values designed for it
Avals = [-2.4724 28.025 9.1404 581.71 3394.8];
Bvals = [1.7978 -6.7881 -197.69 -271.69 -8219.5];
Cvals = [-0.1702 -9.3959 163.53 -319.35 5698.6];
Dvals = [-0.2308 7.1 -64.118 105.98 -1621.9];
Evals = [0.2302 -5.3688 38.124 -46.681 1813.8];
A = polyval(Avals,p_GPa);
B = polyval(Bvals,p_GPa);
C = polyval(Cvals,p_GPa);
D = polyval(Dvals,p_GPa);
E = polyval(Evals,p_GPa);
coefs = [A B C D E];
T_l_K = polyval(coefs,x_mol);