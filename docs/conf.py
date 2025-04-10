#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is the configuration file for the Sphinx documantation build.

@author: thoverga
"""

# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# make shure that the package directory is in the sys path
import os, sys
from pathlib import Path


# Note that on the github workflow the build is executed in the docs folder,
# Thus if that is the case, we need to go up the foldertree
curfolder = os.path.abspath(".")
if "docs" in curfolder:
    # when executing in docs folder
    basefolder = Path(curfolder).parents[0]
else:
    # when executing in basefolder
    basefolder = curfolder

sys.path.insert(0, os.path.join(str(basefolder), "docs"))
sys.path.insert(0, str(basefolder))
sys.path.insert(0, os.path.join(str(basefolder), "src", "valenspy"))

# ValEnsPy must be imported when testing and building the documentation
# locally. However this is overkill for RTD service, so only import it for
# local builds (or on the VSC server)
if "/home/" in str(basefolder) or "/dodrio/" in str(basefolder):
    import valenspy


# logofile = os.path.join(basefolder, "docs", "logo_wide_1280x640.jpeg")


# -- Project information -----------------------------------------------------

project = "ValEnsPy"
copyright = "2024, Cordex-BE team"
# author = "Thomas Vergauwen"


# =============================================================================
# Load extensions list
# =============================================================================

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",  # Autodocument functions
    "sphinx_rtd_theme",  # Use the read the docs theme
    "sphinx.ext.viewcode",  # Button to go to source code
    "sphinx_copybutton",  # Copy button (for examples etc)
    "sphinx.ext.napoleon",  # To convert Numpydocstring to readable format
    "sphinx.ext.autosummary",  # Create neat summary tables
    "myst_parser",  # for including md files (readme)
    "sphinx.ext.autosectionlabel",  # for cross linking
    "nbsphinx",  # to render the notebook examples in the doc
    "sphinx_design",  # for the design of the doc
    "autodoc_diagnostic",  # custom autodoc for the Diagnostic class
    "sphinx.ext.intersphinx",  # Cross-referencing other projects
]

# Configuration for intersphinx to enable cross-referencing
intersphinx_mapping = {
    "xclim": ("https://xclim.readthedocs.io/en/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "datatree": ("https://xarray-datatree.readthedocs.io/en/latest/", None),
    "xesmf": ("https://xesmf.readthedocs.io/en/stable/", None),
}

# =============================================================================
# General configuration
# =============================================================================

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# Turn on sphinx.ext.autosummary
autosummary_generate = True
autosummary_ignore_module_all = True
autosummary_undoc_members = True
numpydoc_class_members_toctree = False
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "private-members": True,
}


# Specify which file formats to render:
source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}
source_encoding = "utf-8"


# The master toctree document.
master_doc = "index"


# When building the doc, sphinx will try to import all the depending packages,
# this is not needed and problematic when building the docs in a clean environment
# So specify which packages can be mocked

autodoc_mock_imports = [
    "xarray",
    "dask",
    "netCDF4",
    "pandas",
    "shapely",
    "xarray-datatree",
    "matplotlib",
    "scipy",
]

# List of patterns, relative to source directory (= docs), that match files and
# directories to ignore when looking for source files.
exclude_patterns = [
    "_build",
    "_templates",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
]

add_function_parentheses = False
add_module_names = False
show_authors = False  # section and module author directives will not be shown
todo_include_todos = False  # Do not show TODOs in docs


# Make sure the target is unique
autosectionlabel_prefix_document = True


# =============================================================================
# Options for HTML output
# =============================================================================
# These options are on the layout of the generated html pages


html_theme = "pydata_sphinx_theme"
html_title = "ValEnsPy documentation"
html_short_title = "ValEnsPy docs"
# html_logo = "logo_small.svg"
html_static_path = ["_static"]
html_css_files = ["custom.css"]  # to specify custom color palets etc
html_show_sphinx = True
html_show_copyright = True
htmlhelp_basename = "ValEnsPy"  # Output file base name for HTML help builder.
html_use_smartypants = True
html_show_sourcelink = True


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/CORDEX-be2/ValEnsPy",
            "icon": "fab fa-github-square fa-xl",
        },
    ],
    # "navbar_center": ["version-switcher", "navbar-nav"],
    # "switcher": {
    #     # The json url must be a full path operationally !!!
    #     "json_url": "https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/docs/_static/custom.css",  # this file contains a dict of all versions to show
    #     "version_match": f"v{version}",  # currently being browsed
    # },
}


# Try to remove the white/dark theme switch (because it is bugging)
# html_context = {
#    "default_mode": "light" #light or dark theme
# }


# html_theme_options["navbar_end"] = ["navbar-icon-links"]

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = f"MetObs Toolkit {metobs_toolkit.__version__} documentation"

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = "MetObs Toolkit documentation"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = "logo_small.svg"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = "_static/logo/favicon.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
# html_css_files = [
#     "custom.css",
# ]
