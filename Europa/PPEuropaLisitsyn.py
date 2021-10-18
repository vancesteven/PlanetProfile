import numpy as np
import config as cfg
bodyname = 'Europa'
bodycode = 502

import spiceypy as spice
spice.furnsh('Utilities/spice/'+cfg.spicePCK)
_, (a,b,c) = spice.bodvcd(bodycode, 'RADII', 3)
R = np.sqrt((a**2 + b**2 + c**2)/3)

### Custom changes made to these profiles

# Overwrites the default values using user-defined lists

# add "field : [profile 1 value, profile 2 value, ...]"
# make sure all lists have the same length
"""
Example 1:
Producing 3 profiles.
1st profile has (Tb_K,ocean_wpct) = (270.0 , 10)
2nd profile has (Tb_K,ocean_wpct) = (271.0 , 5)
3rd profile has (Tb_K,ocean_wpct) = (272.0 , 0)
##############
CustomDict = {
    'Tb_K' : [270.0 , 271.0 , 272.0],
    'ocean_wpct' : [10 , 5 , 0]
}

numProfiles = 3
###############

Example 2:
Using np.linspace and meshgrid to create profiles on a "grid".
We want 25 profiles, one for each pair of
    5 values of Tb_K: [270.0,270.5,271.0,271.5,272.0]
    5 values of ocean_wpct: [0 , 2 , 4 , 6 , 8]
###############
Tb_K_vals = np.linspace(270.0 , 272.0 , 5)
ocean_wpct_vals = np.linspace(0 , 8 , 5)
Tb_K_arr , ocean_wpct_arr = np.meshgrid(Tb_K_vals , ocean_wpct_vals)

CustomDict = {
    'Tb_K' : Tb_K_arr.flatten(),
    'ocean_wpct' : ocean_wpct_arr.flatten()
}

numProfiles = 25
###############
"""

Tb_K_vals = np.linspace(270.0 , 272.0 , 5)
ocean_wpct_vals = np.linspace(0 , 8 , 5)
Tb_K_arr , ocean_wpct_arr = np.meshgrid(Tb_K_vals , ocean_wpct_vals)

CustomDict = {
    'Tb_K' : Tb_K_arr.flatten(),
    'ocean_wpct' : ocean_wpct_arr.flatten()
}

numProfiles = 25

### Construction of default dictionary for Europa

# Orbital and plotting parameters for use in LayeredInduction Response
OrbitalDict = {
    'peaks_Hz' : [4.946e-5, 2.473e-5, 3.259e-6],
    'f_orb' : 2*np.pi/3.55/86400, # frequency of Europa's orbit [rad/s]
    'ionos_bounds' : 100.e3,
    'ionosPedersen_sig' : 30/100e3,
    'ionos_only' : [],
    'PLOT_SIGS' : True,
    'ADD_TRANSITION_BOUNDS' : False
}

# Bulk and surface properties
# note: C/MR2 uses value from Anderson et al. 1998
BulkSurfaceDict = {
    'rho_kgm3' : 2989., # average density [kg/m^3]
    'R_m' : R*1e3, # radius [m]
    'M_kg' : 4.7991e22, # mass [kg]
    'gsurf_ms2' : 1.428,
    'Tsurf_K' : 110,
    'Psurf_MPa' : 0,
    'Cmeasured' : 0.346, # C/MR^2 value
    'Cuncertainty' : 0.005, # C/MR^2 uncertainty
    'FeCore' : True, # boolean for whether planet has an iron core
    'rhoFe' : 8000., # density of pure iron
    'rhoFeS' : 5150. # density of iron sulfate
}

# Mantle heat properties
# cold case
MantleHeatDict = {
    'kr_mantle' : 4., # rock conductivity, (Cammarano et al. 2006, Table 4)
    'Qmantle_Wm2' : 2.2e11 / 4 / np.pi / BulkSurfaceDict['R_m']**2,
    'QHmantle' : 0.,
    'EQUIL_Q' : False
}

# Porosity of the rock
PorosityDict = {
    'POROUS_ROCK' : False,
    'phi_surface' : 0.
}

# Core properties
CoreDict = {
    'xFeS_meteoritic' : 0.0405,
    'xFeS' : 0.55,
    'xFe_core' : 0.0279,
    'XH2O' : 0.0035,
    'rho_sil_withcore_kgm3' : 3644.
}

# Ocean properties (included as sub-dictionary with key 'Ocean')
OceanDict = {
    'ocean_comp' : 'MgSO4', # composition of the ocean
    'ocean_wpct' : 10.0, # % concentration of solute in ocean
    'Tb_K' : 270. # temperature at the bottom of the ice layer
}

# combines all the dictionaries into one containing all the fields
PlanetDict = {
    'name' : bodyname,
    **OrbitalDict,
    **BulkSurfaceDict,
    **MantleHeatDict,
    **PorosityDict,
    **CoreDict,
    **OceanDict
}

### converts a 'dictionary' to a structured array with 'num' repetitions of said dictionary (for multiple profiles)

def dictionaryToStructuredArray(dictionary , num):
    keys = list(dictionary.keys())

    dt = [] # use list instead of numpy since dtype cannot be an array

    for i in range(len(keys)):
        if type(dictionary[keys[i]]) == str:
            dt.append( (keys[i],f'U{len(dictionary[keys[i]])}')) # resolves issue with string datatype in structured array
        else:
            dt.append( (keys[i] , type(dictionary[keys[i]]) ) ) # generate the dtype list for the structured array

    output = np.empty( (num) , dtype=dt ) # preallocate structured array

    for i in range(len(keys)):
        output[0][keys[i]] = dictionary[keys[i]] # create the first profile

    for i in range(1,num):
        output[i] = output[0] # add remaining profiles

    return output

Planet = dictionaryToStructuredArray(PlanetDict , numProfiles)

### amends Planets structured array using custom changes made in first section

for key in CustomDict.keys():
    Planet[key] = CustomDict[key]


### Creating the seismic dictionary
# Attenuation parameters based on those described in Cammarano et al. 2006
IceI = {
    'B_aniso_iceI' : 0.56,
    'gamma_atten_iceI' : 0.2,
    'g_aniso_iceI' : 22.
}
IceII = {
    'B_aniso_iceII' : 0.56,
    'gamma_atten_iceII' : 0.2,
    'g_aniso_iceII' : 25.
}
IceIII = {
    'B_aniso_iceIII' : 0.56,
    'gamma_atten_iceIII' : 0.2,
    'g_aniso_iceIII' : 27.
}
IceV = {
    'B_aniso_iceV' : 0.56,
    'gamma_atten_iceV' : 0.2,
    'g_aniso_iceV' : 28.
}
IceVI = {
    'B_aniso_iceVI' : 0.56,
    'gamma_atten_iceVI' : 0.2,
    'g_aniso_iceVI' : 30.
}
Mantle = {
    'B_aniso_mantle' : 0.56,
    'gamma_atten_mantle' : 0.2,
    'g_aniso_mantle': 30.
}

Seismic = {
    'LOW_ICE_Q' : 1.,
    'QScore' : 1e4,
    'coreEOS' : 'sulfur_core_partition_SE15_1pctSulfur.tab',
    'mantleEOS' : 'CV3hy1wt_678_1.tab',
    **IceI,
    **IceII,
    **IceIII,
    **IceV,
    **IceVI,
    **Mantle
}

### Creating the Params dictionary

Params = {
    'cfg' : cfg,
    'wlims' : [np.log(0.001),np.log(1000)],
    'foursubplots' : True,
    'Legend' : False,
    'LegendPosition' : 'North',
    'ylim' : [910,1230],
    'Pseafloor_MPa' : 350,
    'nsteps_iceI' : 200,
    'nsteps_ocean' : 350,
    'nsteps_ref_rho' : 30,
    'nsteps_mantle' : 500,
    'nsteps_core' : 10,
    'nsteps_colororder' : 'mcbkgrm',
    'Temps' : [250., 252.5, 255., 260., 265., 270., 273.],
    'LineStyle' : '--',
    'wrefLine' : '--',
    'wref' : [0,34]
}