.. _advanced_install:

Advanced installation
=====================

Valenspy is a Python package but has some non-python dependencies. Therefore, the easist way to install ValEnsPy is to use conda.

Dependencies
-------------

ValEnsPy officially supports Python 3.10 and above. 
Its main dependencies are xarray (>v2025.03.0 for DataTree support) and dask. Additionally, ValEnsPy utilizes the existing ecosystem of xarray based weather and climate packages, including xesmf, xclim, and intake-esm. 
Therefore, Valenspy has a non-python dependency through xesmf, namely ESMF (esmpy). However, if regridding functionality is not needed, ValEnsPy can be installed and used without ESMF.

Installing with conda
---------------------

To install with conda ensure that you have either `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ or `Anaconda <https://docs.continuum.io/free/anaconda/>`__ installed, then run the following command in your terminal:

.. code-block:: shell

    #WIP still needs to be published to conda-forge

This will install ValEnsPy and all its dependencies, including ESMF (esmpy) if it is not already installed in the environment.

.. _install.pip:

Installing with pip
-------------------

Valenspy can be installed via pip from `PyPI <https://pypi.org/project/ValEnsPy/>`__. 

.. warning::
    Installing ValEnsPy with pip will not install ESMF (esmpy). If you require regridding functionality, either install esmpy in your environment seperately or use conda to install ValEnsPy.

    .. code-block:: shell

        conda install -c conda-forge esmpy

In your terminal run the following command:

.. code-block:: shell
    
    pip install ValEnsPy

Installing on Windows
---------------------

ValEnsPy is not yet fully supported on Windows because the xesmf package (in particular the ESMF dependency) is not supported on Windows.
If you are using Windows and do not require regridding functionality, you can :ref:`install ValEnsPy using pip <_install.pip>` without installing ESMF (esmpy) separately. If the regridding functionality is we recommend using a Linux or MacOS environment.

Installation from source
------------------------

#WIP - link to developer guide/contributing pages

Testing the installation
------------------------

#WIP




