* This directory contains:
  1. Python code for running PlanetProfile
  1. *Default* settings to use for each body (in the Default/Body/PPBody.py files)
  1. *Default* general settings.
  1. Templates for user setting files to be copied the the top-level directory.

Files in this directory should be changed only for development purposes. Do `python -m PlanetProfile.install PPinstall` in the top-level directory to copy default profiles to working directories for experimentation. This also copies the config files to the top-level directory.

*All config files in the top-level directory override settings in the default files here. All PlanetProfile/PlanetProfile/Default/Body/PPBody.py files will be overridden by the files found at PlanetProfile/Body/PPBody.py.*
