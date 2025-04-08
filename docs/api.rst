.. _api_docs:

.. TODO - add the pre-made diagnostics to the toctree

=============
API reference
=============

This page gives an overview of all ValEnsPy objects, functions, and methods. All classes and functions exposed in :ref:`valenspy`.* are public API.

Modules
-----------

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:

   input
   processing
   diagnostic

Pre-made Diagnostic
-------------------
Most users will use pre-made diagnostics. They are organized into four categories and listed below.

Model2Self
==========

See :py:class:`~diagnostic.Model2Self` for all shared functionality.

.. automodule:: diagnostic._model2self
   :members:
   :undoc-members:

Model2Ref
=========

See :py:class:`~diagnostic.Model2Ref` for all shared functionality.

.. automodule:: diagnostic._model2ref
   :members:
   :undoc-members:
   :imported-members:

Ensemble2Self
=============

See :py:class:`~diagnostic.Ensemble2Self` for all shared functionality.

.. automodule:: diagnostic._ensemble2self
   :members:
   :undoc-members:

Ensemble2Ref
============

See :py:class:`~diagnostic.Ensemble2Self` for all shared functionality.

.. automodule:: diagnostic._ensemble2ref
   :members:
   :undoc-members:
