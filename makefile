# For installing and running PlanetProfile.
# USAGE:
#	make: Print available command line arguments.
#	make install: Copy necessary files for running PlanetProfile.
#	make uninstall: Remove files outside the PlanetProfile directory that were placed by 'install'.
#	make clean: Remove files ignored by GitHub (output data and figures).
#	make pp: Open Matlab with default arguments.
#	make <body>: Same as pp, but opens and runs the input file for the named body.
#
# DESCRIPTION:
#	Various commands for setting up and running PlanetProfile. Designed for accessibility for new users and ease of use for advanced users.
# 
# AUTHOR: Marshall Styczinski (mjstyczi@uw.edu), 2018-07-03

SHELL := /bin/bash

refprop=0

mbodies="Callisto" "Enceladus" "Europa" "Ganymede" "Titan"
figs="figures"
ppdir=$(shell pwd)

uname=$(shell uname -s)
ifeq ($(uname),Darwin)
	# Mac OS
	matlabpath=$(shell matlabp=($$(find /Applications/MATLAB*/toolbox/local -type d)) ; echo $$matlabp)
else
	# Other Unix
	matlabpath=$(shell matlabp=($$(find $$HOME/usr/local/MATLAB/*/toolbox/local -type d)) ; echo $$matlabp)
endif

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
	@echo "pp:		Open Matlab with default arguments."
	@echo "<body>:		Same as pp, but opens and runs the input file"
	@echo "		  for the named body."
	@echo 

clean:
	@bodies=($(mbodies)) ; \
	echo "Clearing .gitignore files for:" $${bodies[@]}  ; \
	for body in $${bodies[@]} ; do \
		flist1=$$(cat $$body/.gitignore) ; \
		flist2=$$(cat $$body/$(figs)/.gitignore) ; \
		for line in $$flist1 ; do \
			rm -f $$body/$$line ; \
		done ; \
		for line in $$flist2 ; do \
			rm -f $$body/$(figs)/$$line ; \
		done ; \
	done

pp:
	@echo "WIP: Not finished yet."

install:
ifeq ($(refprop),1)
	mkdir -p /opt/refprop
	cp Thermodynamics/librefprop.so-master/librefprop.dylib /opt/
	mkdir -p /opt/refprop/fluids /opt/refprop/mixtures
	cp Thermodynamics/librefprop.so-master/files/*.fld /opt/refprop/fluids/
	cp Thermodynamics/librefprop.so-master/files/*.mix /opt/refprop/mixtures/
endif

ifneq ($(uname),Darwin)
	@echo "Warning: If your Matlab install location differs from the default"
	@echo "  at ~/usr/local/MATLAB, you will need to edit the makefile"
	@echo "  to match."
endif

	@pathdirs=($$(find * -type d -not -path *version* -not -path *input*)) ; \
	if [ -z $(matlabpath)/startup.m ] ; then \
		echo "cd $(ppdir)" > startup.m ; \
		for subdir in $${pathdirs[@]} ; do \
			echo "addpath('$$subdir')" >> startup.m ; \
		done ; \
		mv startup.m $(matlabpath)/ ; \
		echo "Matlab startup file created with PlanetProfile folders in path:" ; \
		echo "  $(matlabpath)/startup.m" ; \
	else \
		echo "cd $(ppdir)" >> $(matlabpath)/startup.m ; \
		for subdir in $${pathdirs[@]} ; do \
			echo "addpath('$$subdir')" >> $(matlabpath)/startup.m ; \
		done ; \
		echo "PlanetProfile folders added to Matlab path in startup file:" ; \
		echo "  $(matlabpath)/startup.m" ; \
	fi
	@echo " "
	@echo "Installation complete. In Matlab, run startup command, or close and"
	@echo "  relaunch, to automatically add PlanetProfile dirs to your Matlab path."

uninstall:
ifeq ($(refprop),1)
	rm /opt/librefprop.dylib
	rm /opt/refprop/fluids/*.fld
	rm /opt/refprop/mixtures/*.mix
	rmdir /opt/refprop/fluids /opt/refprop/mixtures	
	rmdir /opt/refprop
	@echo "Refprop files removed."
endif

	@pathdirs=($$(find * -type d -not -path *version* -not -path *input*)) ; \
	echo "cd $(ppdir)" > startup.m ; \
	for subdir in $${pathdirs[@]} ; do \
		echo "addpath('$$subdir')" >> startup.m ; \
	done
	@if diff startup.m $(matlabpath)/startup.m >/dev/null 2>&1 ; then \
		rm $(matlabpath)/startup.m ; \
		echo "Matlab startup.m removed." ; \
	else \
		echo "$(matlabpath)/startup.m contains more than just the PlanetProfile path. Matlab startup.m not modified." ; \
	fi
	rm startup.m
	@echo " "
	@echo "Uninstall complete."
	@echo "Delete this directory and all subdirectories to finish purge."
	@echo " "
