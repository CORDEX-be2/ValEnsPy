:html_theme.sidebar_secondary.remove:


ValEnsPy
============================

A Python package to validate gridded climate model data.




.. toctree::
   :hidden:
   :maxdepth: 7

   Home <self>
   Getting started <getting_started>
   Examples <doc_examples/index>
   API <API_doc>
   Contributing <contribution_pages/index>

.. container:: button

    :doc:`Getting started <getting_started>` :doc:`Contributing <contribution_pages/index>`
    :doc:`Examples <doc_examples/index>`



Description
-----------

The package is designed to validate gridded climate model data.
The building-blocks are xarray and dask, in particular xarray-DataTrees are used to manage multiple models and gridded observations.
It contains three main components:
* **Input processor**: Reads the input data and prepares it into a CF-compliant format.
* **Pre-processor**: Transforms all input datasets to enable intercomparison.
* **Diagnostic**: Calculates the diagnostic metrics and plots the results.

The package is designed to be modular, allowing users to use the components independently or together.

A visual representation of the package structure is shown below.

.. image:: package_structure.png
    :alt: logo
    :width: 700


Indices and tables
----------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
