# Installation

Currently, the package is in development and not yet available on [PyPI](https://pypi.org/).
Therefore, to install the package locally, follow the instructions below.
For usage on the VSC (Tier 1) platform, see the [VSC installation guide](INSTALL_VSC.md).

## Local installation

### Required software

The following software (or equivalent) is required to set up a developer environment:
* [Anaconda](https://anaconda.org/)

> [!TIP]
> Miniconda is small and easy to [install](https://docs.anaconda.com/free/miniconda/#quick-command-line-install).

Make sure you have this software installed before proceeding:
```bash
which conda #This should return the path to the conda executable
```

### Setup a developer environment
1. Clone the ValEnsPy repository locally:

> [!TIP]
> First time using git/github? Read this small git [guide](GIT.md) to get started.

```bash
git clone git@github.com:CORDEX-be2/ValEnsPy.git
```

2. Create a conda environment and install the required packages.

You can specify the location that you want to install the conda environment. The default location is in the home directory, which is often limited in size.
```bash
conda config --add envs_dirs /path/to/your/envs
```

Then create the environment with python and poetry.
```bash
conda create -n valenspy_dev python=3.11 esmpy poetry=1.8 -c conda-forge
source activate valenspy_dev
```
Now install the required packages using poetry: All required packages are listed in the pyproject.toml file and will be installed automatically by poetry.

> [!WARNING]
> Make sure you are in the branch you want to work in before installing the packages, as different branches may have different dependencies. The 'dev' branch is a good starting point.

```bash
cd ValEnsPy
git checkout dev
```

Finally install all required packages.

> [!NOTE]
> Optionally useful packages can be installed by installing specific groups of packages by adding the `--with` flag followed by the group name.
> Available groups are:
> - `examples` to be able to run example notebooks (e.g. for ipykernel for jupyter notebooks)
> - `docs` for building the documentation
> - `dev` for development tools

```bash
#Recommended for users
poetry install --with examples
#For developers
#poetry install --with examples --with dev --with docs
```

To test if the installation was successful, run the tests:

```bash
python tests/push_tests/import_test.py
```

3. Create a branch for local development which is a copy of the **dev** branch.

If this using git and branching is new to you, please read the [git guide](GIT.md) first.

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

> [!TIP]
> Some IDEs (like VSCode) have a built-in git interface that can be used to commit and push changes.

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

## Additional information
For additional developer guidelines see the Notes for developers page in the Contributing section of the documentation.
