# Contributing VSC

All contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

## Ways to contribute
There are many ways to contribute to the ValEnsPy package:
1. Partake in one of the [discussions](https://github.com/CORDEX-be2/ValEnsPy/discussions) or create a new one.
2. Create a new [issue](https://github.com/CORDEX-be2/ValEnsPy/issues) to help us improve the package.
  - New feature requests
  - Bug Reports

3. Contribute to the codebase. This can be done by:
  - Taking on an open issue
  - Adding a new diagnostic (see [Adding a new diagnostic](#adding-a-new-diagnostic)

4. Help improve the documentation. This can be done by:
  - Adding a new example
  - Improving the general documentation

## Contribute to the codebase
Ready to make code contributions? Here's how to set up ValEnsPy for local development.
Here you can find a guide for setting up a local development environment on the **VSC (Vlaams Supercomputer Centrum) specifically**.

### Required software

The following software (or equivalent) is required to set up a developer environment:
* [Anaconda](https://anaconda.org/)

On the VSC Miniconda can be used instead of Anaconda.
To install Miniconda on the VSC, follow the instructions on the [VSC website](https://docs.vscentrum.be/software/python_package_management.html#install-miniconda).

To test if Miniconda is installed correctly, the following command should provide the path to the conda executable ($VSC_data/miniconda3/bin/conda):
```bash
which conda
```

### Using the existing conda environment
On the VSC, there is an existing conda environment that can be used for development. This environment is called `valenspy_dev` and can be activated using the following command:

TODO Test if this works on the VSC!?

### Setup a developer environment

0. Choose a location for your development environment. Ideaaly this should be within the 2022_200 project directory.
For example:

```bash
cd /dodrio/scratch/projects/2022_200/project_output/__INSTITUTE__/__VSC_USERNAME__
```

1. Clone the ValEnsPy repository locally:

```bash
git clone git@github.com:CORDEX-be2/ValEnsPy.git
```
TODO: Add link to Git guide for cloning the repository (using SSH keys).

2. Create a conda environment and install the required packages.

First initialize the conda environment and install python (version 3.9) and poetry (a python package manager).

```bash
conda create -p ./valenspy_dev python==3.9 poetry
source activate valenspy_dev
```
Note that the conda environment is created in the current directory (./valenspy_dev) - this is optional and can be changed to a different location.

Now install the required packages using poetry: All required packages are listed in the pyproject.toml file and automatically installed by poetry.

```bash
cd ValEnsPy
git checkout dev #Or the branch you want to work on (e.g. structure)
poetry install
```

3. Create a branch for local development which is a copy of the **dev** branch.:

```bash
# checkout the dev branch
git checkout dev
git pull

# Create a new local branch and switch to it.
git branch name-of-your-bugfix-or-feature
git checkout name-of-your-bugfix-or-feature
```
Now you can make local changes.

5. Push your code online:

```bash
# Add your changes to your commit
git add -A
# Write commit text
git commit -m "Some text describing your code changes in this commit"
# Push your branch online
#only the first time:
git push --set-upstream origin name-of-your-bugfix-or-feature
#all other times
git push
```

## Pull Request Guidelines
Once your branch has been *pushed* to github, you can create a *Pull request* in github. Make sure that you have **referred the corresponding issues** to the *Pull request*.
If your code adaptations are still *work-in-progress* add the ![Static Badge](https://img.shields.io/badge/WIP%20-%20%23A21079) label to it. For each push, github will perform a list of checks (package building, version control, functionality test, os-tests, documentation build test), in order to merge your contributions these tests must all be successful.

If your code is ready for review, you can add the ![Static Badge](https://img.shields.io/badge/Ready_for_Review%20-%20%230315E4) label to it.

After the code review, and all review marks are resolved, your contributions will be merged to the *dev* branch.

## Versioning/Tagging
From time to time the *dev* branch will be merged with the master with a new [*Release tag*](https://github.com/CORDEX-be2/ValEnsPy/releases).

## Acknowledgement
This file is inspired by the [RavenPy](https://github.com/CSHS-CWRA/RavenPy) project. Thank you for the inspiration!”.
