.. _getting_started:

Getting started
=================

A guide to enable you to start productively using ValEnsPy as efficiently as possible.

Installation
------------

Valenspy is a pure python package but it has some non-python dependencies (such as ESMpy). The easiest way to install ValEnsPy is to use the conda package manager. If you prefer using pip ensure that you have the required dependencies installed. For help with the install check out the :ref:`advanced installation <advanced_install>` page.

.. grid:: 1 2 2 2
    :gutter: 4

    .. grid-item-card:: Working with conda?
        :class-card: install-card
        :columns: 12 12 6 6
        :padding: 3

        Using Conda

        ++++++++++++++++++++++

        .. code-block:: bash

            #TODO: Add conda install command for non developers

    .. grid-item-card:: Prefer pip?
        :class-card: install-card
        :columns: 12 12 6 6
        :padding: 3

        Using pip

        ++++

        .. code-block:: bash

            #Ensure non-python dependencies (ESMpy) are installed, e.g. using conda
            #conda install -c conda-forge esmpy
            pip install valenspy

    .. grid-item-card:: In-depth instructions?
        :class-card: install-card
        :columns: 12
        :padding: 3

        Installing a specific version? Installing from source? Check the advanced
        installation page.

        +++

        .. button-ref:: advanced_install
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            Learn more


Why ValEnsPy?
-------------

WIP - 2 sentence description of the package and its purpose.

Still not convinced? Check out the :ref:`Why ValEnsPy <why_ValEnsPy>` page.

What do I need to get started?
------------------------------

WIP - 3 sentences about knowledge needed to get started with ValEnsPy (xarray!)

Got what it takes? Check out the :ref:`quick overview <quick-overview>` page to learn the key concepts of ValEnsPy or check out the :ref:`examples <examples_index>`.

.. Add a list of examples here similar to pandas.


.. toctree::
   :maxdepth: 2
   :hidden:

   Getting started <self>
   why_ValEnsPy
   advanced_install
   quick-overview
   faq       