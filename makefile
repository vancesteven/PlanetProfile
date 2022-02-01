# For installing and running PlanetProfile.
# USAGE:
#	make: Print available command line arguments.
#	make install: Copy necessary files for running PlanetProfile.
#	make uninstall: Remove files outside the PlanetProfile directory that were placed by 'install'.
#	make clean: Remove files ignored by GitHub (output data and figures).
#
# DESCRIPTION:
#	Various commands for setting up and running PlanetProfile.
#	Designed for accessibility for new users and ease of use for advanced users.
# 
# AUTHOR: Marshall J. Styczinski (itsmoosh@gmail.com), 2018-07-03
#	Last updated 2022-01-01

SHELL := /bin/bash

# Set to 0 for PlanetProfile versions using Refprop
refprop=1

mbodies="Ariel" "Callisto" "Enceladus" "Europa" "Ganymede" "Miranda" "Titan" "Triton" "Test"
figs="figures"
cdpp=cd $(shell pwd)

foundmatlab=1
windows=1

# Check for likely install locations
ifeq ($(shell [ -d /Applications/MATLAB* ] ; echo $$?),0)
	# Mac OS
	matlabdir=/Applications/MATLAB*
	foundmatlab=0
else ifeq ($(shell [ -d /mnt/c/Program\ Files/MATLAB* ] ; echo $$?),0)
	# Windows
	windows=0
	matlabdir=/mnt/c/Program\ Files/MATLAB/*
	driveltr=$(shell ppdir=$$(pwd) ; echo $${ppdir:5:1} | tr "[:lower:]" "[:upper:]")
	rempath=$(shell ppdir=$$(pwd) ; echo $${ppdir:6})
	cdpp=cd $(driveltr):$(rempath)
	foundmatlab=0
else ifeq ($(shell [ -d /mnt/c/Program\ Files\ \(x86\)/MATLAB* ] ; echo $$?),0)
	# Windows
	windows=0
	matlabdir=/mnt/c/Program\ Files\ \(x86\)/MATLAB/*
	driveltr=$(shell ppdir=$$(pwd) ; echo $${ppdir:5:1} | tr "[:lower:]" "[:upper:]")
	rempath=$(shell ppdir=$$(pwd) ; echo $${ppdir:6})
	cdpp=cd $(driveltr):$(rempath)
	foundmatlab=0
else ifeq ($(shell [ -d /usr/local/MATLAB/* ] ; echo $$?),0)
	# Linux
	matlabdir=/usr/local/MATLAB/*
	foundmatlab=0
else ifeq ($(shell [ -d $$HOME/usr/local/MATLAB/* ] ; echo $$?),0)
	# Linux
	matlabdir=$$HOME/usr/local/MATLAB/*
	foundmatlab=0
else
	# Non-standard installation location. Make the user find the right spot.
	foundmatlab=1
endif
matlabsupath=$(shell find $(matlabdir)/toolbox/local -maxdepth 0 -type d)

default:
	@echo "Pass arguments to 'make' for one of the functions described below."
	@echo "Your command line argument should look like: make command"
	@echo 
	@echo "<no command>:	Print this list."
	@echo "install:	Copy necessary files for running PlanetProfile"
	@echo "		  and set default Matlab path."
	@echo "uninstall:	Remove files outside the PlanetProfile directory"
	@echo "		  that were placed by 'install'."
	@echo "clean:		Remove files ignored by GitHub (output data and figures)."
	@echo 

clean:
	@bodies=($(mbodies)) ; \
	echo "Clearing data and figure files for:" $${bodies[@]} ; \
	for body in $${bodies[@]} ; do \
		rm -f $$body/*.dat $$body/*.txt $$body/*.mat $$body/*.asv ; \
		rm -f $$body/$(figs)/*.eps $$body/$(figs)/*.fig $$body/$(figs)/*.png $$body/$(figs)/*.pdf ; \
	done

install:
# Copy Refprop files and libraries into system folders
ifeq ($(refprop),0)
	mkdir -p /opt/refprop
	cp Thermodynamics/librefprop.so-master/librefprop.dylib /opt/
	mkdir -p /opt/refprop/fluids /opt/refprop/mixtures
	cp Thermodynamics/librefprop.so-master/files/*.fld /opt/refprop/fluids/
	cp Thermodynamics/librefprop.so-master/files/*.mix /opt/refprop/mixtures/
endif

	@# Get list of PlanetProfile directories containing important files
	@#   and print them into a file named startup.m
	@# First, check if SeaFreeze is installed and grab it if not
	@SeaFinstalled=$$(find Thermodynamics/SeaFreeze/* | wc -l) ; \
	if [[ $$SeaFinstalled -eq 2 ]] ; then \
		echo "SeaFreeze has not been installed yet." ; \
		read -p "Would you like to install it now using git? ([Y]/n)" doSeaFinstall ; \
		if [[ "$$doSeaFinstall" =~ ^([nN][oO]|[nN])$$ ]] ; then \
			SeaFinstalled=-1 ; \
		else \
			git clone git@github.com:Bjournaux/SeaFreeze.git Thermodynamics/SeaFreeze/SeaFreeze && SeaFinstalled=0 ; \
		fi ; \
	else \
		SeaFinstalled=0 ; \
	fi ; \
	pathdirs=($$(find * -type d -not -path *version* -not -path *.git* -not -path *Python*)) ; \
	echo "$(cdpp)" > startup.m ; \
	for subdir in $${pathdirs[@]} ; do \
		echo "addpath('$$subdir')" >> startup.m ; \
	done ; \
	echo "addpath(fullfile('$(shell pwd)/Utilities/spice/mice/src/mice/'))" >> startup.m
	echo "addpath(fullfile('$(shell pwd)/Utilities/spice/mice/lib/'))" >> startup.m
	if [[ $$SeaFinstalled -eq -1 ]] ; then \
		echo "WARNING: SeaFreeze was not installed when the startup path" ; \
		echo "  file was created. Delete startup.m from the Matlab" ; \
		echo "  toolbox/local directory and run make install again after" ; \
		echo "  installing SeaFreeze in order to add it to the startup path." ; \
	fi

ifeq ($(foundmatlab),1)
	@echo "WARNING: Your Matlab install location was not found."
	@echo "Copy the file named startup.m from the PlanetProfile"
	@echo "  directory to the Matlab toolbox/local directory"
	@echo "  to complete installation."
	@echo
	@echo "Next, in Matlab, go to:"
	@echo "  Preferences -> Matlab -> General -> Toolbox path caching"
	@echo "  and click on the button for Update Toolbox Path Cache."
	@echo "Then, at the Matlab command prompt, run the startup command, or close and"
	@echo "  relaunch, to automatically add PlanetProfile dirs to your Matlab path."
else ifeq ($(windows),0)
	@# Windows folders are most often not writable. Make the user do it.
	@echo "To complete installation, copy startup.m from the PlanetProfile"
	@echo "  directory into your Matlab toolbox/local directory at:"
	@echo "  $(matlabsupath)"
	@echo " "
	@echo "Next, in Matlab, go to:"
	@echo "  Preferences -> Matlab -> General -> Toolbox path caching"
	@echo "  and click on the button for Update Toolbox Path Cache."
	@echo "Then, at the Matlab command prompt, run the startup command, or close and"
	@echo "  relaunch, to automatically add PlanetProfile dirs to your Matlab path."
else
	@echo "Matlab toolbox/local found: $(matlabsupath)"

	@# Check if the Matlab startup file exists, and if so, if it's just for PlanetProfile
	@#   so we can avoid duplicates.
	@if cmp -s startup.m "$(matlabsupath)"/startup.m ; then \
		echo "Matlab startup file for PlanetProfile is already in place:" ; \
		echo "  $(matlabsupath)/startup.m" ; \
	elif [ ! -f "$(matlabsupath)"/startup.m ] ; then \
		cp startup.m "$(matlabsupath)"/startup.m ; \
		echo "Attempted to create Matlab startup file with PlanetProfile dirs in path:" ; \
		echo "  $(matlabsupath)/startup.m" ; \
	else \
		cat startup.m >> "$(matlabsupath)"/startup.m ; \
		echo "PlanetProfile dirs appended to Matlab path in startup file:" ; \
		echo "  $(matlabsupath)/startup.m" ; \
	fi
	@rm startup.m
	@echo " "
	@echo "Installation finished. If errors ocurred above or if startup.m is"
	@echo "  still in the PlanetProfile directory, paste into the Matlab"
	@echo "  toolbox/local/ dir."
	@echo " "
	@echo "Next, in Matlab, go to:"
	@echo "  Preferences -> Matlab -> General -> Toolbox path caching"
	@echo "  and click on the button for Update Toolbox Path Cache."
	@echo "Then, at the Matlab command prompt, run the startup command, or close and"
	@echo "  relaunch, to automatically add PlanetProfile dirs to your Matlab path."
endif

uninstall:
# Delete Refprop files we placed
ifeq ($(refprop),0)
	rm /opt/librefprop.dylib
	rm /opt/refprop/fluids/*.fld
	rm /opt/refprop/mixtures/*.mix
	rmdir /opt/refprop/fluids /opt/refprop/mixtures	
	rmdir /opt/refprop
	@echo "Refprop files removed."
endif

ifeq ($(foundmatlab),1)
	@echo "WARNING: Your Matlab install location was not found."
	@echo "Delete the PlanetProfile lines added to your Matlab"
	@echo "  startup file at MATLAB/toolbox/local/startup.m"
	@echo "  and delete this repository to complete uninstall."
else
	@# Recreate startup.m file to check against the one Matlab is using
	@pathdirs=($$(find * -type d -not -path *version* -not -path *input*)) ; \
	echo "$(cdpp)" > startup.m ; \
	for subdir in $${pathdirs[@]} ; do \
		echo "addpath('$$subdir')" >> startup.m ; \
	done
	
	@# If the Matlab startup.m file matches the one we just made, it exists only
	@#   for PlanetProfile. Delete it to return the system to its initial state.
	@# Otherwise, startup.m has been modified by the user, so we should not touch it.
	@if cmp -s startup.m "$(matlabsupath)"/startup.m ; then \
		rm "$(matlabsupath)"/startup.m ; \
		echo "Attempted to remove Matlab startup.m." ; \
	elif [ ! -f "$(matlabsupath)"/startup.m ] ; then \
		echo "There was no Matlab startup file in $(matlabsupath)." ; \
		echo "Was PlanetProfile installed before uninstalling?" ; \
	else \
		echo "$(matlabsupath)/startup.m contains more than just the PlanetProfile path." ; \
		echo "Matlab startup.m not modified." ; \
	fi
	
	@# Get rid of our checking file
	@rm startup.m
	
	@echo
	@echo "Uninstall complete, assuming no errors occurred above."
	@echo "Delete this directory and all subdirectories to complete uninstall."
	@echo
endif
