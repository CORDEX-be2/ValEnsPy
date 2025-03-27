:html_theme.sidebar_secondary.remove:


ValEnsPy
============================

ValEnsPy makes working with gridded climate and weather data easier, reproducible and collaborativeby covering the entire workflow from data processing to diagnostics and from single model evaluations up to multiple ensemble comparisons. 

**Useful links**:
`Code Repository <https://github.com/CORDEX-be2/ValEnsPy>`__ |
`Issues <https://github.com/CORDEX-be2/ValEnsPy/issues>`__ |
`Releases <https://github.com/CORDEX-be2/ValEnsPy/releases>`__

.. grid:: 1 2 2 2
    :gutter: 4
    :padding: 2 2 0 0
    :class-container: sd-text-center

    .. grid-item-card:: Getting started
        :img-top: _static/index_getting_started.svg
        :class-card: intro-card
        :shadow: md

        Eager to use ValEnsPy? Check out the getting started guides.

        +++

        .. button-ref:: getting_started
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the getting started guides.

    .. grid-item-card::  User guide
        :img-top: _static/index_user_guide.svg
        :class-card: intro-card
        :shadow: md

        The user guide provides in-depth information on the
        key concepts of pandas with useful background information and explanation.

        +++

        .. button-ref:: 
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the user guide

    .. grid-item-card::  API reference
        :img-top: _static/index_api.svg
        :class-card: intro-card
        :shadow: md

        The reference guide contains a detailed description of
        the pandas API. The reference describes how the methods work and which parameters can
        be used. It assumes that you have an understanding of the key concepts.

        +++

        .. button-ref:: API_doc
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the reference guide      

    .. grid-item-card::  Developer guide
        :img-top: _static/index_contribute.svg
        :class-card: intro-card
        :shadow: md

        Saw a typo in the documentation? Want to improve
        existing functionalities? The contributing guidelines will guide
        you through the process of improving pandas.

        +++

        .. button-ref:: contribution_pages/index
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            To the development guide

.. toctree::
   :hidden:
   :maxdepth: 7

   Home <self>
   Getting started <getting_started>
   Examples <doc_examples/index>
   API <API_doc>
   Contributing <contribution_pages/index>

.. Description
.. -----------

.. The package is designed to validate gridded climate model data.
.. The building-blocks are xarray and dask, in particular xarray-DataTrees are used to manage multiple models and gridded observations.
.. It contains three main components:
.. * **Input processor**: Reads the input data and prepares it into a CF-compliant format.
.. * **Pre-processor**: Transforms all input datasets to enable intercomparison.
.. * **Diagnostic**: Calculates the diagnostic metrics and plots the results.

.. The package is designed to be modular, allowing users to use the components independently or together.

.. A visual representation of the package structure is shown below.

.. .. image:: package_structure.png
..     :alt: logo
..     :width: 700


.. Indices and tables
.. ----------------------

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
