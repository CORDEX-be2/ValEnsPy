Notes for developers
======================

Here is a page dedicated for developers to share notes and instructions for
other developers.

On general, the `deployment` folder contains executable (bash) scripts, to
streamline some typical development procedures.

Merging and Pull Requests
---------------------------

Once you are satisfied with your code and are ready to share it with the rest of the team, you can create a *Pull request* in github.
A *Pull request* is a request to merge your code into the *dev* branch. It allows other developers to review your code and suggest changes.
Even if you are not ready to merge your code, you can still create a *Pull request*. This is useful to have an overview of incoming changes.

Here are some guidelines for creating a *Pull request*:
- Make sure that you have referred the corresponding issues to the *Pull request* (if any).
- If your code adaptations are still *work-in-progress*, add the ![Static Badge](https://img.shields.io/badge/WIP%20-%20%23A21079) label to it - this will indicate to other developers that the code is not ready for review.
- If your code is ready for review, you can add the ![Static Badge](https://img.shields.io/badge/Ready_for_Review%20-%20%230315E4) label to it.

After the code review, and all review marks are resolved, your contributions will be merged to the *dev* branch.

Versioning/Tagging
---------------------

From time to time the *dev* branch will be merged with the master with a new [*Release tag*](https://github.com/CORDEX-be2/ValEnsPy/releases).

Packaging
------------

ValEnsPy is desinged to be packaged. This is done automatically on github
by a workflow, if:

* Code is pushed to `main`
* Code is pushed to `dev`
* 'test' is found in a commit message

You can also build the package, (especially handy for testing), on your local machine.
To do this, you run the `dev_pipeline.sh` script:

.. code-block::

    source deployment/dev_pipeline.sh



Handling Dependencies
-------------------------

If you make code adaptations and introduce a new dependency (aka import a new package),
you need to specify this in order for packaging. This is to make shure that all
dependencies are campatible with each other.

If you introduce a new dependency, you must add it in the `deployment/dev_pipeline.sh` script
under "Add dependencies to pyproject".


.. code-block::

    #example: add shapely as a dependency, add this in the dev_pipeline script
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

The `docs/` folder contains all the file to build the documentation. To view the
documentation, go to `docs/_build` and open the `index.html` file in your browser.
(If the `docs/_build` is empty, you need to build the documentation, see below)


The documentation is a combination of

* Static file: These are standard files, that are converted to HTML. Most of them are RST and MD files,
and if you do not change them, the documentation will not change as well (trivial).

* Dynamic part: When you make a code contribution, your contribution must be represented in the API page of the documentation.
So this means that each time you make a code contribution, the documentation (only the API part), must be upated as well. This is the
reason that every developer must know how to build the documentation form there working branch.

To build the documentation, simply run the `deployment/dev_pipeline.sh` file.
If you want to build the documentation without building the package, you must
install all the dependencies for the documentation (see `deployment/dev_pipeline.sh`),
and build it using sphinx.

.. code-block::

    cd docs
    sphinx-build -a -E -v . ./_build/

or execute the `docs/sphinx_build` script.


When building the documentation, take a look at the sphinx errors and warnings,
they indicate syntax and formatting issues.


Code-Standards
----------------------

@kobe: You should descide on these standards


* Docstring: Use Numpy-style docstrings, other types are not rendered well in the documentation. See https://numpydoc.readthedocs.io/en/latest/example.html#example




Dev-checklist
---------------

* Adding a new (non-standard) dependencies

  #. Add the package in the `deployment/dev_pipeline.sh` (see Handling Dependencies)

  #. run the `deployment/dev_pipeline.sh`

  #. Check if the dependency is written in the `pyproject.toml` file. (Do not make changes in it.)


* Adding a new module (adding new files to automatically document).

  #. Add the module in the toctree of `docs/API_doc`

* Pre-commit

  #. Run Pre-commit before merging for code formatting and linting.
