Notes for developers
======================

Here is a page dedicated for developers to share notes and instructions for
other developers.

On general, the `deployment` folder contains executable (bash) scripts, to
streamline some typical development procedures.


Packaging
------------

ValEnsPy is desinged to be packaged. This is done automatically on github
by a workflow, if:

* Code is pushed to `main`
* Code is pushed to `dev`
* 'test' is found in a commit message

You can also build the package, (especially handy for testing), on your local machine.
To do this, you run the `build_package.sh` script:

.. code-block::

    source deployment/build_package.sh



Handling Dependencies
-------------------------

If you make code adaptations and intoduce a new dependency (aka import a new package),
you need to specify this in order for packaging. This is to make shure that all
dependencies are campatible with each other.

If you introduce a new dependency, you must add it in the `deployment/build_package.sh` script
under "Add dependencies to pyproject".


.. code-block::

    #example: add shapely as a dependency, add this in the build_package script
    poetry add shapely


In general it is best to not specify the version of the depency manually, poetry
takes care of this.


If you added a dependency, that is not required for using the ValEnsPy package
(i.g. Poetry -- for packaging, jupyter-lab -- for making examples, ...), you can
add the dependencies in a seperated group. This is to make shure that if a user
whants to install the package, these dependencies are not installed. But a developer
has the option to install them by

.. code-block::

    #install the package (and the required dependencies)
    poetry install

    #install the package with extra groups
    poetry install --with dev, packaging

    # install the package with all dependencies
    poetry install --all-extras






Documentation builds
--------------------------

bla ba



CF-Convention
---------------

bla bla





