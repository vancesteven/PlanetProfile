# PlanetProfile, Software framework for constructing interior structure models based on planetary properties. 
# Copyright (C) 2018  Steven D. Vance
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

import os, sys, numpy as np
from PlanetProfile.Main import run

if __name__ == '__main__':
    # Command line args
    nArgs = len(sys.argv)
    clArg = None
    fNames = None
    if nArgs > 1 and ('PP' not in sys.argv[1] and '.txt' not in sys.argv[1]):
        # Body name was passed as command line argument
        bodyname = sys.argv[1]

        # Additional command line arguments
        if nArgs > 2:
            if 'PP' in sys.argv[2]:
                print('PP in CL arg 2 -- interpreting as (list of) filename(s).')
                fNames = sys.argv[2:]
            elif '.txt' in sys.argv[2]:
                print('.txt in CL arg 2 -- interpreting as (list of) filename(s) to reload.')
                fNames = sys.argv[2:]
                clArg = 'reload'
            else:
                clArg = sys.argv[2]
                if nArgs > 3:
                    if 'PP' in sys.argv[3]:
                        fNames = sys.argv[3:]
                    else:
                        print(f'Too many command line args passed. Ignoring command "{sys.argv[3]}" and any after it.')
    elif 'PP' in sys.argv[1]:
        print('PP in first CL arg -- interpreting as (list of) filename(s).')
        fNames = sys.argv[1:]
        if len(fNames) == 1:
            bodyname = os.path.split(fNames[0])[0]
            fNames = os.path.split(fNames[0])[1:]
        else:
            bodynames = [os.path.split(fName)[0] for fName in fNames]
            if np.all(bodynames == bodynames[0]):
                bodyname = bodynames[0]
                fNames = [os.path.split(fName)[1] for fName in fNames]
            else:
                bodyname = ''
    elif '.txt' in sys.argv[1]:
        print('.txt in first CL arg -- interpreting as (list of) filename(s) to reload.')
        fNames = sys.argv[1:]
        bodyname = ''
        clArg = 'reload'
    else:
        # No command line argument, ask user which body to run
        bodyname = input('Please input body name: ')

    run(bodyname=bodyname, opt=clArg, fNames=fNames)
