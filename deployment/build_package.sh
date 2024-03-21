#!/bin/bash

#This script will make a local build of the package on a Unix system.


# ----- Setup paths ---------
cd ..
WORKDIR=$(pwd)
DOCDIR=${WORKDIR}/docs
DISTDIR=${WORKDIR}/dist
SRCDIR=${WORKDIR}/src/valenspy


# ---- Add dependencies to pyproject  ------

#add packages that 
poetry add xarray
poetry add pandas

# --- Build the package ---------

#remove previous builds
cd ${DISTDIR}
rm *.whl
rm *.tar.gz

cd ${WORKDIR} #(maybe this is not needed)
poetry update #to update the poetry.lock file (aka use the latest, valid, dependencies)
poerty install --all-extras
poetry show #print out some dep. information

cd ${DISTDIR}
poetry build #build the package 




# ---- Use the build to run the testscripts ---------



