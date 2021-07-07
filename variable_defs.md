Description
===========

This file lists the variable names and definitions used in the `PPBody.py` files.

Planet Object
=============

The `Planet` object is a [numpy structured array](https://numpy.org/doc/stable/user/basics.rec.html), containing fields describing the body. The specific fields used depend on the properties of the body.

There are `numProfiles` values for each field, used by `PlanetProfile.py` to construct models for multiple possible setups at once.

Orbital and plotting parameters
-------------------------------
- `peaks_Hz`  (list of floats):
- `f_orb` (float): the frequency of the body's orbit \[1/s\]
- `wlims` (list of floats):
- `ionos_bounds` (float):
- `ionosPedersen_sig` (float):
- `ionos_only` (??):
- `PLOTS_SIGS` (boolean):
- `ADD_TRANSITION_BOUNDS` (boolean):

Bulk and surface properties
---------------------------
- `rho_kgm3` (float): average density of the body \[kg/m$^3$\]  
- `R_m` (float): radius of the body \[m\]  
- `M_kg` (float): mass of the body \[kg\]  
- `gsurf_ms2` (float):  
- `Tsurf_K` (float):  
- `Psurf_MPa` (float):  
- `Cmeasured` (float): $C/MR^2$ (polar moment of inertia of body, normalized to $MR^2$)
- `Cuncertainty` (float): the uncertainty in the measurement of $C/MR^2$

Mantle Heat Properties
----------------------
- `kr_mantle` (float??): rock conductivity
- `Qmantle_Wm2` (float):
- `QHmantle` (??):
- `EQUIL_Q` (??)

Core Properties
---------------
- `FeCore` (boolean): True iff the body has an iron core
- `rhoFe` (float): density of pure iron
- `rhoFeS` (float): density of iron sulfate
- `xFeS_meteoritic` (float):
- `xFeS` (float): proportion of Iron Sulfide in the core (??)
- `xFe_core` (float): proportion of pure iron in the core (??)
- `XH2O` (float): proportion of water in the body (where?) (??)
- `rho_sil_withcore_kgm3` (float): density of silicon (??)
- `phi_surface` (??):

Ocean Properties
----------------
- `ocean_comp` (str): solute found in ocean
- `ocean_wpct` (float): % weight of solute found in ocean
- `Tb_K` (float): temperature at bottom of ice shell


Seismic Object
==============

The `Seismic` object is a [python dictionary](https://docs.python.org/3/tutorial/datastructures.html#dictionaries), containing fields describing seismic properties of the body and relevant materials. The specific fields used depend on the properties of the body.

Attenuation Parameters
--------------------------
Replace `N` with a roman numeral representing the ice phase $(N \in \{I,II,III,V,VI\})$:
- `B_aniso_iceN` (float): 
- `gamma_atten_iceN` (float):
- `g_aniso_iceN` (float):

Similar parameters exist for the mantle:
- `B_aniso_mantle` (float): 
- `gamma_atten_mantle` (float):
- `g_aniso_mantle` (float):

Equations of State
------------------
- `coreEOS` (str): name of file containing equation of state for the core
  - only used if the body has an iron core
- `mantleEOS` (str): name of file containing equation of state for the mantle

Other (??)
----------
- `LOW_ICE_Q` (??):
- `Qscore` (??):


Params Object
=============
The `Partams` object is a [python dictionary](https://docs.python.org/3/tutorial/datastructures.html#dictionaries), containing fields describing a variety of miscellaneous properties of the model, including information regarding plots. The specific fields used depend on the properties of the body.

- `cfg` (object): a config object containing fields with variety of relevant 'configuration' information
- `wlims` (list of floats): (?? seems to be the same as Planet.wlims with orbital parameters)
- `foursubplots` (boolean?):
- `Legend` (boolean): decides whether plots will have legends?
- `LegendPosition` (str): decides the location of the legend in the plots
- `ylim` (2-length list of floats):
- `Pseafloor_MPa` (float): pressure at the seafloor \[MPa\]
- `Temps` (list of floats): (??)
- `LineStyle` (str): style of lines in plot?
- `wrefLine` (str):
- `wref` (list of ??): 

nsteps fields
-------------
Represent the number of steps used by `PlanetProfile.py` in variety of layers
- `nsteps_iceI` (int): ice shell
- `nsteps_ocean` (int): ocean
- `nsteps_ref_rho` (int): ??
- `nsteps_mantle` (int): mantle
- `nsteps_core` (int): core
- `nsteps_colororder` (str): used to plot the colors of all the layers