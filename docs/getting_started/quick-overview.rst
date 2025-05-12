.. _quick-overview:

Quick overview
==============

Here we illustrate the key concepts used in ValEnsPy. 

Components
----------

Valenspy consists of three main components:
- **Input**: Gathering raw data, loading it and transforming it to ValEnsPy complaint xarray DataSet or DataTree with uniform naming conventions.
- **Processing**: User dependent processing steps, e.g. regridding, masking, etc.
- **Diagnostics**: The computation and visualization of the diagnostics.

The main data structures used in ValEnsPy are xarray DataSets and DataTrees. In particular, the new xarray DataTree structure is used to manage (multiple) ensembles of gridded models and observations.
The image below illustrates the main components of ValEnsPy and how they interact with each other.

.. image:: /_static/images/valenspy_overview.png
    :alt: Overview diagram of ValEnsPy
    :align: center
    :width: 80%

Input
^^^^^

The input component is responsible for gathering the raw data, loading it and transforming it to ValEnsPy complaint xarray DataSet or DataTree with uniform naming conventions.

The latter is done by the `InputConvertors` class, which is essentially a dataset specific pre-processing function applied when loading the data. It essentially consists of:

- A dictionary to map the raw data variables to the `CMIP6-CORDEX <https://docs.google.com/spreadsheets/d/1qUauozwXkq7r1g-L4ALMIkCNINIhhCPx/edit?gid=1672965248#gid=1672965248>`_ variable names.
- If required, a dataset specific function to transform the raw data to the ValEnsPy compliant format.

For standard datasets, ValEnsPy has built in input processors but users can also easily define their own input processors.

Gathering and loading the data is done by the `Manager` class, which creates a catalog of all available datasets and utilizes the `intake-esm` to make that data searchable and loadable.
The creation of the catalog is semi-automatic, i.e. only the base directory and a pattern for the files need to be specified. On shared shared machines this could only have to be done once after which the catalog can be used by all users. The catalog is stored in a `yaml` file and can be easily shared with others. The catalog is then used to load the data into xarray DataSets or DataTrees.
When loading the data, the `Manager` class also applies a set of pre-processing steps to each respective dataset through the aformentioned `InputConvertors`.

Processing
^^^^^^^^^^

The processing component enables users to apply the required processing steps to the data. These are a combination of simple xarray operations time selection, masking and more complex operations like regridding and calculating indicators.
Where required ValEnsPy extends existing processing functionality in particular to support the new xarray DataTree structure.

Diagnostics
^^^^^^^^^^^

Finally, the diagnostics are used to compute and visualize the results. Each diagnostic represents a diagnsotic function and at least one plot.
Diagnostics are categorized into 4 groups, each with slightly different scope and functionality:

- **Single model diagnostics (Model2Self)**: Compute and visualize aspects of a single model. e.g: spatial average of a variable.
- **Single model to reference diagnostics (Model2Ref)**: Compute and visualize aspects of a single model compared to a reference dataset. e.g: is the average spatial bias of a variable.
- **Multi-model diagnostics (Model2Model)**: These diagnostics compute and visualize aspects of an whole ensemble of models. e.g: The "most extreme" model in the ensemble with respect to some variables.
- **Multi-model to reference diagnostics (Model2Ref)**: These diagnostics compute and visualize aspects of an whole ensemble of models compared to a reference dataset. e.g: The spread of the bias of the ensemble.

The diagnostics functions are applied on the xarray DataSets or DataTrees resulting in some form of output (pandas DataFrame, xarray DataSet or DataTree, dictionary, etc.) which can be saved or visualized with the diagnostic plot functions.

Within ValEnsPy there are `some prexisting diagnostics <_api_docs>`_, but users can also define their own diagnostics.