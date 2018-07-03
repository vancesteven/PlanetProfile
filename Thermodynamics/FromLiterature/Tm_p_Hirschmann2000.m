function Tm_oC = Tm_p_Hirschmann2000(P_GPa)
coeffs = [-5.140,132.899,1120.661];
Tm_oC = polyval(coeffs,P_GPa);