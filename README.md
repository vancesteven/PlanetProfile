# PlanetProfile
Matlab software for constructing 1D interior structure models based on planetary properties. Self-consistent thermodynamics are used for fluid, rock, and mineral phases. Sound speeds, attenuation, and electrical conductivities are computed as outputs.

The main code is called from an input file containing all the planetary data.  Ideally, no tweaks to the main code are needed in order to change the outputs of the model.  

TODOs:
Modularization is not complete. 
The mantle temperature profile is currently fixed to be the cold case from Cammarano et al. 2006
MgSO4 is the only currently available equation of state
NH3-H2O is in progress
NaCl and Na2SO4 have data sets and preliminary equations of state, but melting behavior has not been coded 
