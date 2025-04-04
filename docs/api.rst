.. _api_docs:

=============
API reference
=============

This page gives an overview of all ValEnsPy objects, functions, and methods. All classes and functions exposed in `valenpsy.*` are public. 
As the package is still in its infancy, the API may change in the future.

Top-level functions
===================

.. Remove remap_cdo in the future

.. autosummary::
   :toctree: _autosummary

   processing.remap_xesmf
   processing.remap_cdo
   processing.convert_units_to
   processing.xclim_indicator

.. Diagnostics
.. ===========

.. Constructors
.. ------------

.. .. autosummary::
..    :toctree: _autosummary

..    Diagnostic
..    Model2Self
..    Model2Ref
..    Ensemble2Self
..    Ensemble2Ref

.. Attributes
.. ----------

.. .. autosummary::
..    :toctree: _autosummary

..    Diagnostic
.. ..   Diagnostic.description
.. ..   Diagnostic.diagnostic_function
.. ..   Diagnostic.plotting_function
..    .. Model2Self.plot_type
..    .. Ensemble2Ref.iterative_plotting

.. Functions
.. ---------

.. .. autosummary::
..    :toctree: _autosummary

..    Diagnostic
.. ..   Diagnostic.apply
.. ..   Diagnostic.plot
.. ..   Model2Self.apply_dt

.. .. Unify so that only apply is used (not apply_dt)

Pre-made diagnostics
====================

Model2Self
----------

.. automodule:: diagnostic._model2self
   :toctree: _autosummary
   :members:
   :undoc-members:
   :imported-members:
   :template: custom-module-template.rst

Model2Ref
---------

.. automodule:: diagnostic._model2ref
   :toctree: _autosummary
   :members:
   :undoc-members:
   :imported-members:
   :template: custom-module-template.rst

Ensemble2Self
-------------

.. automodule:: diagnostic._ensemble2self
   :toctree: _autosummary
   :members:
   :undoc-members:
   :imported-members:
   :template: custom-module-template.rst

Ensemble2Ref
------------

.. automodule:: diagnostic._ensemble2ref
   :toctree: _autosummary
   :members:
   :undoc-members:
   :imported-members:
   :template: custom-module-template.rst

