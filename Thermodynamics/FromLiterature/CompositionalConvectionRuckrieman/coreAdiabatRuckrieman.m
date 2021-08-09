function Tcore_K = coreAdiabatRuckrieman(nTbs,nsteps_core,Tcmb,r_core_m)
% from Ruckrieman et al. 2018 (Equation 1)
alpha_c = 9e-5; % core thermal expansivity Ruckrieman2018 Table 1
c_pc = 800; % core heat capacity Ruckrieman2018 Table 1
G = 6.67408e-11; % gravitational constant
rho_c_avg = 5000; % average density in core
Tcore_K = zeros(nTbs,nsteps_core);
for iT = 1:nTbs
    Tcore_K(iT,:) = Tcmb(iT)*exp(-alpha_c/c_pc*2/3*pi*G*rho_c_avg*(r_core_m(iT,:).^2-r_core_m(iT,1)^2));
end