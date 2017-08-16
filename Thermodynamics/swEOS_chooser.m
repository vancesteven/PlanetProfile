function swEOS = swEOS_chooser(EOStype)
% swEOS    Seawater Equation of state
%=========================================================================
%
% USAGE:  swEOS = swEOS_chooser(EOStype)
% DESCRIPTION:
%    Set the seawater equation of state to be used by the column-convect
%    model.  Input units for all functions:
%    Temperature: Kelvin
%    Pressure: Bar, absolute
%    Salinity: molal for MgSO4, absolute salinity (g/kg) for GSW
%
% INPUT: 
%   EOStype (string): mgso4 | gsw
%
% OUTPUT:
%   handles to relevant functions
% 
% AUTHOR:  Steve Vance (svance@jpl.nasa.gov)
%
% REFERENCES:  
% Marion et al (2005), Geochimica et Cosmochimica Acta, Vol. 69, No. 2, pp. 259?274
% Vance and Brown (2005), doi:10.1016/j.icarus.2005.06.005
% Vance and Brown (2013), Geochimica et Cosmochimica Acta, 110:176?189
% Vance et al. (2014), Planetary And Space Science, 96:62?70
% McDougall and P. M. Barker (2011) SCOR/IAPSO WG, 127:1?28
%=========================================================================

gsw302path = '' ;
gsw305path = '' ;
% CALLER: general purpose
% CALLEE:  mgso4_loader polyvaln
switch(EOStype)
    case('mgso4')
        swEOS.adtg=@mgso4_adtg; % (S,T,P) [K/bar] 
        swEOS.alpha=@mgso4_alpha; % (S,T,P) [1/K] 
        swEOS.beta=@mgso4_beta; % (S,T,P) [1/mol]
        swEOS.C=@Larionovp01molalMgSO4; % (T,P) [S/m]
        swEOS.cp=@mgso4_cp; % (S,T,P) [J/kg-K]
        swEOS.dens=@mgso4_dens; % (S,T,P) [kg/m^3]
        swEOS.tfreezing=@mgso4_t_freezing; % (S,P) [K]
        swEOS.latent_heat_melting=@mgso4_latent_heat_melting; % L = mgso4_latent_heat_melting(S,P)
        swEOS.pden=@mgso4_pden; % (S,T,P,PR) [kg/m^3] 
        swEOS.ptmp=@mgso4_ptmp; % (S,T,P,PR) [K}
        swEOS.vel=@mgso4_vel; % (S,T,P) [km/s] 
    case('gsw302')
        rmpath gsw305path
        addpath gsw302path
        % mod functions convert T(K) and P(bar) inputs to t(oC) and p(dbar)
        % inputs
        swEOS.adtg=@modgsw_adiabatic_lapse_rate_t_exact; % (SA,t,p) [K/bar]
        swEOS.alpha=@modgsw_alpha_wrt_t_exact; % (SA,t,p) [1/K]
        swEOS.beta=@modgsw_beta_const_t_exact; % (SA,t,p) [ kg/g ]
        swEOS.C = @modgsw_C_from_SA; % (SA,t,p) [S/m]
        swEOS.cp=@modgsw_cp_t_exact; % (SA,t,p) [J/(kg*K)]
        swEOS.dens=@modgsw_rho_t_exact; % (SA,t,p) [ kg/m^3 ]
        swEOS.latent_heat_melting=@modgsw_latentheat_melting; % (SA,p) [J/kg]
        swEOS.pden=@modgsw_pot_rho_t_exact; % (SA,t,p,p_ref) [kg/m^3]
        swEOS.ptmp=@modgsw_pt_from_t; % (SA,t,p,p_ref) [oC]
        swEOS.tfreezing=@modgsw_t_freezing; % (SA,p) [oC]; modgsw_t_freezing below returns T in K instead 
        swEOS.pfreezing=@modgsw_pressure_freezing_CT; % (SA,p) [bar]; modgsw_t_freezing below returns T in K instead 
        swEOS.vel=@modgsw_sound_speed_t_exact; % (SA,t,p) [m/s]
        swEOS.vel_ice = @mgsw_sound_speed_ice; % (t,p) [m/s]
        swEOS.velP_iceIh = @modgsw_sound_speed_ice; % (t,p) [m/s]
        swEOS.velS_iceIh = @velS_gsw; % (t,p) [m/s]
    case('gsw305')
%         rmpath gsw305path
%         addpath gsw302path
        %below is what was there before. check to see if the above breaks
        %things.
%         rmpath ../gsw_matlab_v3_02 ../gsw_matlab_v3_02/library
%         addpath ../gsw_matlab_v3_05_5 ../gsw_matlab_v3_05_5/library ../gsw_matlab_v3_05_5/thermodynamics_from_t
        % mod functions convert T(K) and P(bar) inputs to t(oC) and p(dbar)
        % inputs
        swEOS.adtg=@modgsw_adiabatic_lapse_rate_from_t; % (SA,t,p) [K/bar]
        swEOS.alpha=@modgsw_alpha_wrt_t_exact; % (SA,t,p) [1/K]
        swEOS.beta=@modgsw_beta_const_t_exact; % (SA,t,p) [ kg/g ]
        swEOS.C = @modgsw_C_from_SA; % (SA,t,p) [S/m]
        swEOS.cp=@modgsw_cp_t_exact; % (SA,t,p) [J/(kg*K)]
        swEOS.dens=@modgsw_rho_t_exact; % (SA,t,p) [ kg/m^3 ]
        swEOS.latent_heat_melting=@modgsw_latentheat_melting; % (SA,p) [J/kg]
        swEOS.pden=@modgsw_pot_rho_t_exact; % (SA,t,p,p_ref) [kg/m^3]
        swEOS.ptmp=@modgsw_pt_from_t; % (SA,t,p,p_ref) [oC]
        swEOS.tfreezing=@modgsw_t_freezing; % (SA,p) [oC]; modgsw_t_freezing below returns T in K instead 
        swEOS.pfreezing=@modgsw_pressure_freezing_CT; % (SA,p) [bar]; modgsw_t_freezing below returns T in K instead 
        swEOS.vel=@modgsw_sound_speed_t_exact; % (SA,t,p) [m/s]
        swEOS.velP_iceIh = @modgsw_sound_speed_ice; % (t,p) [m/s]
        swEOS.velS_iceIh = @velS_gsw; % (t,p) [m/s]
end

% Convert T from K to oC, P from absolute bars to atmospheric gauge
% decibars for use with gsw.
function adtg = modgsw_adiabatic_lapse_rate_t_exact(SA,T,P)
% for compatibility with gsw_v3_02
adtg = gsw_adiabatic_lapse_rate_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5))*1e5;  
% Returns value in K/bar for compatibility with mgso4.

function adtg = modgsw_adiabatic_lapse_rate_from_t(SA,T,P)
adtg = gsw_adiabatic_lapse_rate_from_t(SA,T-gsw_T0,10*(P-gsw_P0/1e5))*1e5;  
% Returns value in K/bar for compatibility with mgso4.

function alpha = modgsw_alpha_wrt_t_exact(SA,T,P)
alpha = gsw_alpha_wrt_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5));

function beta = modgsw_beta_const_t_exact(SA,T,P)
beta = gsw_beta_const_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5));

function C = modgsw_C_from_SA(SA,T,P)
SP = gsw_SP_from_SA(SA,10.1325,0,0); % set lat lon to 0
if SP<2 | SP>42
    error('SP out of bounds')
end
C = 0.1*gsw_C_from_SP(SP,T-gsw_T0,10*(P-gsw_P0/1e5));

function cp = modgsw_cp_t_exact(SA,T,P)
cp = gsw_cp_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5));

function dens = modgsw_rho_t_exact(SA,T,P)
dens = gsw_rho_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5));

function latent_heat_melting = modgsw_latentheat_melting(SA,P)
latent_heat_melting = gsw_latentheat_melting(SA,10*(P-gsw_P0/1e5));

function pden = modgsw_pot_rho_t_exact(SA,T,P,P_ref)
pden = gsw_pot_rho_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5),10*(P_ref-gsw_P0/1e5));

function ptmp = modgsw_pt_from_t(SA,T,P,P_ref)
ptmp_oC = gsw_pt_from_t(SA,T-gsw_T0,10*(P-gsw_P0/1e5),10*(P_ref-gsw_P0/1e5));
ptmp = ptmp_oC + gsw_T0;

function P = modgsw_pressure_freezing_CT(SA,T)
P = fzero(@(P) T-gsw_T0-gsw_t_freezing(SA,P),[0 10000]);
% P = fzero(@(P) gsw_chem_potential_water_t_exact(SA,T-gsw_T0,P)-...
%     gsw_chem_potential_water_ice(T-gsw_T0,P),[0 21000]); % minimize gibbs energy over the 210 MPa range of the Feistel 2006 EOS.
P = 0.1*P+gsw_P0/1e5; % pressure in bar

function T = modgsw_t_freezing(SA,P)
% t_oC = fzero(@(T) gsw_chem_potential_water_t_exact(SA,T,10*(P-gsw_P0/1e5))-...
%     gsw_chem_potential_water_ice(T,10*(P-gsw_P0/1e5)),[-100 20]); % minimize gibbs energy over the 210 MPa range of the Feistel 2006 EOS.
t_oC = gsw_t_freezing(SA,10*(P-gsw_P0/1e5)); % t_freezing = gsw_t_freezing(SA,p,saturation_fraction)
T = t_oC+gsw_T0;

function vel = modgsw_sound_speed_t_exact(SA,T,P)
vel = gsw_sound_speed_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5));

function vel = modgsw_sound_speed_ice(T,P)
vel = gsw_sound_speed_ice(T-gsw_T0,10*(P-gsw_P0/1e5));

function vel = velS_gsw(T,P)
vP = gsw_sound_speed_ice(T-gsw_T0,10*(P-gsw_P0/1e5)); %m/s
KS = gsw_kappa_ice(T-gsw_T0,10*(P-gsw_P0/1e5)); % 1/Pa
V = gsw_specvol_ice(T-gsw_T0,10*(P-gsw_P0/1e5)); % m3/kg
vel = sqrt(0.75*(vP.^2-V./KS)); % m/s
