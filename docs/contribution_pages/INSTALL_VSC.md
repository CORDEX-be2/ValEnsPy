# Installation

Currently, the package is in development and not yet available on [PyPI](https://pypi.org/).
Therefore, to install the package locally, follow the instructions below.
This quide is specifically for the VSC. 

## Local installation

### Required software

The following software (or equivalent) is required to set up a developer environment:
* [Anaconda](https://anaconda.org/)

We recommend using Miniconda, to install Miniconda follow the instructions on the [VSC website](https://docs.vscentrum.be/software/python_package_management.html#install-miniconda).

> [!TIP]
> Their documentation recommends installing conda on $VSC_DATA but this is very slow, therefore, do the install on the project scratch.

Make sure you have this software installed before proceeding:
```bash
which conda #This should return the path to the conda executable
```

### Using the existing conda environment
On the VSC, there is an existing conda environment that can be used. This environment is called `valenspy_dev`.
This environment can be added to your conda environment list by adding the following path to your .conda/environments.txt file:

```bash
echo "/dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/valenspy_dev_shared" >> ~/.conda/environments.txt
```

You can now activate the environment:
```bash
source activate /dodrio/scratch/projects/2022_200/project_output/RMIB-UGent/valenspy_dev_shared
```

This is ideaal for testing the package, and developing new features if no additional packages are required.
If additional packages are required, a new conda environment should be created - see the next section.

### Setup a developer environment

0. Choose a location for your development environment. Ideaaly this should be within the 2022_200 project directory.
For example:

```bash
cd /dodrio/scratch/projects/2022_200/project_output/__INSTITUTE__/__VSC_USERNAME__
```

1. Clone the ValEnsPy repository locally:

> [!TIP]
> First time using git/github? Read this small git [guide](docs/contribution_pages/GIT.md) to get started.

```bash
git clone git@github.com:CORDEX-be2/ValEnsPy.git
```

2. Create a conda environment and install the required packages.

First initialize the conda environment and install python (version 3.9) and poetry (a python package manager).

```bash
conda create -p ./valenspy_dev python==3.9 poetry
source activate valenspy_dev
```
> [!NOTE]
> The conda environment is created in the current directory (./valenspy_dev) - this is optional and can be changed to a different location.

Now install the required packages using poetry: All required packages are listed in the pyproject.toml file and will be installed automatically by poetry.

> [!WARNING]
> Make sure you are in the branch you want to work in before installing the packages as different branches may have different dependencies. The 'dev' branch is a good starting point.

```bash
cd ValEnsPy
git checkout dev
poetry install --all-extras
```

To test if the installation was successful, run the tests:

```bash
python tests/push_tests/import_test.py
```

3. Create a branch for local development which is a copy of the **dev** branch.

If this using git and branching is new to you, please read the [git guide](docs/contribution_pages/GIT.md) first.

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