#!/bin/bash

#This script will make a local build of the package on a Unix system.
# 1 Add dependencies to the pyproject file
# 2 remove previous builds
# 3 update dependencies and build package
# 4 build the documentation

# ----- Setup paths ---------
#Run this from the root of the project
WORKDIR=$(pwd)
DOCDIR=${WORKDIR}/docs
DISTDIR=${WORKDIR}/dist
SRCDIR=${WORKDIR}/src/valenspy
TESTDIR=${WORKDIR}/tests


# ---- Add dependencies to pyproject  ------

# !!!!!!!!!!!!!!!!!!!!!!!!!
# @ALL: If you introduce a new dependency, add it to this list as
#   poetry add PYPI_NAME_OF_DEPENDENCY
# !!!!!!!!!!!!!!!!!!!!!!!!!

poetry add xarray
poetry add pandas
poetry add shapely
poetry add xarray-datatree
poetry add matplotlib
poetry add cdo
poetry add cf-xarray
poetry add regionmask
poetry add geopandas

#add to dev (= development) group
poetry add --group dev 'pre-commit'


#add to the packaging group
poetry add --group packaging 'poetry'

#add to docs (= documentation) group
poetry add --group docs 'Sphinx'
poetry add --group docs 'sphinx_rtd_theme'
poetry add --group docs 'sphinx_copybutton'
poetry add --group docs 'myst_parser'
poetry add --group docs 'nbsphinx'
poetry add --group docs 'pydata-sphinx-theme'

#add to the examples group
poetry add --group examples gcsfs
poetry add --group examples zarr
poetry add --group examples jinja2

# ==============================================================================
# Build the package
# ==============================================================================

#remove previous builds
mkdir -p ${DISTDIR}
cd ${DISTDIR}
rm *.whl
rm *.tar.gz

cd ${WORKDIR} #(maybe this is not needed)
poetry update #to update the poetry.lock file (aka use the latest, valid, dependencies)
poetry install --with dev --with docs
poetry show #print out some dep. information

cd ${DISTDIR}
poetry build #build the package
cd ${WORKDIR}

# ==============================================================================
# Build the documentation
# ==============================================================================
cd ${DOCDIR}
source sphinx_build
cd ${WORKDIR}

# ==============================================================================
# Run tests
# ==============================================================================

# ---- Use the build to run the testscripts ---------

cd ${TESTDIR}
poetry run python push_tests/import_test.py
cd ${WORKDIR}
