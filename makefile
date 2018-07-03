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

mbodies = "Callisto" "Enceladus" "Europa" "Ganymede" "Titan"
figs = "figures"

uname=$(shell uname -s)
ifeq ($(uname),Darwin)
	# Mac OS
	matlabpath=/Applications/MATLAB*/bin
else
	# Other Unix
	matlabpath=$$HOME/usr/local/MATLAB/*/bin
endif

default:
	@echo "Pass arguments to 'make' for one of the functions described below."
	@echo "Your command line argument should look like: make command"
	@echo 
	@echo "<no command>:	Print this list."
	@echo "install:	Copy necessary files for running PlanetProfile."
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
	@if [ -z $$(which matlab) ] ; then \
		echo "matlab command not found. Run make installpp and be sure the matlab executable is found in your \$$PATH." ; \
		stop ; \
	fi
	matlab -r "folder='.'"
	@echo "WIP: Not finished yet."

install:
	@if [ -z $$(which matlab) ] ; then \
		echo " " >> $$HOME/.bash_profile ; \
		echo "# Added by PlanetProfile" >> $$HOME/.bash_profile ; \
		echo "export PATH=\$$PATH:$$(echo $(matlabpath))" >> $$HOME/.bash_profile ; \
	fi

	mkdir -p /opt/refprop
	cp Thermodynamics/librefprop.so-master/librefprop.dylib /opt/
	mkdir -p /opt/refprop/fluids /opt/refprop/mixtures
	cp Thermodynamics/librefprop.so-master/files/*.fld /opt/refprop/fluids/
	cp Thermodynamics/librefprop.so-master/files/*.mix /opt/refprop/mixtures/
	@echo " "
	@echo "To complete installation, relaunch Terminal or type the following command:"
	@echo "	source ~/.bash_profile"
	@echo " "

uninstall:
	rm /opt/librefprop.dylib
	rm /opt/refprop/fluids/*.fld
	rm /opt/refprop/mixtures/*.mix
	rmdir /opt/refprop/fluids /opt/refprop/mixtures	
	rmdir /opt/refprop
	@echo " "
	@echo "Uninstall complete. Files within this directory have not been modified."
	@echo "Delete this directory and all subdirectories to finish purge."
	@echo "You may also want to delete the lines inserted into your ~/.bash_profile."
	@echo " "