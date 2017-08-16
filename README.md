# PlanetProfile
Matlab software for constructing 1D interior structure models based on planetary properties. Self-consistent thermodynamics are used for fluid, rock, and mineral phases. Sound speeds, attenuation, and electrical conductivities are computed as outputs.

The main code is called from an input file containing all the planetary data.  Ideally, no tweaks to the main code are needed in order to change the outputs of the model.  

Some calculations use Matlab's Parallel Computing package.  If you don't have access to this package then parfor loops should be changed to for loops.  A future version will check and do this automatically.

Calculations with seawater solutions use the Gibbs Seawater package for Matlab: http://www.teos-10.org/software.htm#1

Calculations with NH3 solutions use REFPROP and require a compiled dynamic library based on the REFPROP source code.  The source can be obtained from the National Institute of Standards and Technology https://www.nist.gov/refprop
Access to REFPROP functions is through python 3 using librefprop.so: https://github.com/jowr/librefprop.so
The python capabilities of are employed using the included matlab code refproppy.m
REFPROP version 10, expected in October 2017, will provide matlab functions and Mac modules, which may eliminate the need for the above workarounds.

Rock properties are from Perple_X: http://www.perplex.ethz.ch/
Input files were developed by Fabio Cammarano. Version 6.7.8 is currently being used.

TODOs:
Modularization is not complete. 
Further equations of state are under development
Update to work with REFPROP V10
