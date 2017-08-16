# PlanetProfile
Matlab software for constructing 1D interior structure models based on planetary properties. Self-consistent thermodynamics are used for fluid, rock, and mineral phases. Sound speeds, attenuation, and electrical conductivities are computed as outputs.

The main code is called from an input file containing all the planetary data.  Ideally, no tweaks to the main code are needed in order to change the outputs of the model.  

Some calculations use Matlab's Parallel Computing package.  If you don't have access to this package then parfor loops should be changed to for loops.  A future version will check and do this automatically.

Calculations with NH3 solutions use REFPROP and require a compiled dynamic library based on the REFPROP source code.  The source can be obtained from the National Institute of Standards and Technology https://www.nist.gov/refprop
Access to REFPROP functions is through python 3 using librefprop.so: https://github.com/jowr/librefprop.so
The python capabilities of are employed using the included matlab code refproppy.m
REFPROP version 10, expected in October 2017, will provide matlab functions and Mac modules, which may eliminate the need for the above workarounds.

TODOs:
Modularization is not complete. 
The mantle temperature profile is currently fixed to be the cold case from Cammarano et al. 2006.
MgSO4 is the only currently available equation of state.
NH3-H2O is in progress.
NaCl and Na2SO4 have data sets and preliminary equations of state, but melting behavior has not been coded.
