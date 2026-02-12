# ValEnsPy
A Python package to validate ensembles gridded model data.

## Cloning the repository

To clone the repository, run the following command in your terminal:

```bash
git clone git@github.com:CORDEX-be2/ValEnsPy.git
```

## Installation

Install the development version of the package using conda:

Create a conda environment with the required non python packages
```bash
conda create -n valenspy_dev python=3.11 esmpy poetry=1.8 -c conda-forge
source activate valenspy_dev
```

Install the required packages using poetry. 
> [!IMPORTANT] 
> Ensure that you are in the ValEnsPy directory and on the dev branch (or your branch of choice).
```bash
poetry install --with examples
```

For more detailed installation instructions - [General install guide](docs/contribution_pages/INSTALL.md)

Interested in contributing? Take a look at our [contributing guidelines](docs/contribution_pages/CONTRIBUTING.md).

## Documentation
Take a look at our [documentation](https://cordex-be2.github.io/ValEnsPy/)

## Package structure

The package is structured as follows

![Package structure](docs/package_structure.png) 
