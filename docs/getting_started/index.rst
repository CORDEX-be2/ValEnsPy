.. _getting_started:

Getting started
=================

A guide to enable you to start productively using ValEnsPy as efficiently as possible.

.. _install:

Installation
------------

ValEnsPy has ESMF as a non-python dependency. Therefore, when installing ValEnsPy with pip, ensure that ESMF (esmpy) is already installed in the environment if you wish to use regridding functionality. 

.. warning::
    To install ValEnsPy on Windows see the :ref:`advanced installation page <advanced_install>`.

.. grid:: 1 2 2 2
    :gutter: 4

    .. grid-item-card:: Working with conda?
        :class-card: install-card
        :columns: 12 12 6 6
        :padding: 3

        Using Conda (Recommended)

        ++++++++++++++++++++++

        .. code-block:: bash

            #conda install command for non developers
            #WIP

    .. grid-item-card:: Prefer pip?
        :class-card: install-card
        :columns: 12 12 6 6
        :padding: 3

        Using pip

        ++++

        .. code-block:: bash

            pip install valenspy

    .. grid-item-card:: In-depth instructions?
        :class-card: install-card
        :columns: 12
        :padding: 3

        Installing on Windows? Installing from source or with pip? Check the advanced
        installation page.

        +++

        .. button-ref:: advanced_install
            :ref-type: ref
            :click-parent:
            :color: secondary
            :expand:

            Advanced installation


Why ValEnsPy?
-------------

By utilizing the exisisting xarray ecosystem, ValEnsPy provides a flexible and powerful framework for working with gridded climate and weather data from data processing to diagnostics and from single model evaluations up to multiple ensemble comparisons.

Still not convinced? Check out the :ref:`Why ValEnsPy <why_ValEnsPy>` page.

What do I need to get started?
------------------------------

- A working `installation <advanced_install>`_ of ValEnsPy.
- A basic understanding of `xarray <https://docs.xarray.dev/en/stable/getting-started-guide/index.html>`_ in particular the newly introduced `DataTree <https://docs.xarray.dev/en/stable/user-guide/data-structures.html#datatree>`_ functionality.
- A basic understanding of `pandas <https://pandas.pydata.org/docs/getting_started/index.html>`_

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