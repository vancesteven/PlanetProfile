% Reproduces the dashed line from Ruckrieman 2018 Figure 10a 

p_GPa = 6;
x_eut = 0.11 + 0.187*exp(-0.065*p_GPa);
x_s_below = 0:x_eut/100:x_eut;

T_K_below = IronRichLiquidusRuckrieman(p_GPa,x_s_below);

x_s_above = x_eut:(0.3636-x_eut)/100:0.3636;

T_K_above = SulfurRichLiquidusRuckrieman(p_GPa,x_s_above);

plot( [x_s_below x_s_above] , [T_K_below T_K_above] )
xlabel( "Sulfur concentration (wt%)" )
ylabel( "Liquidus Temperature (K)" )