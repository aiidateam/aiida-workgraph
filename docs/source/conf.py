# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import shutil
from pathlib import Path

sys.path.insert(0, os.path.abspath("../.."))


# -- Project information -----------------------------------------------------

project = "AiiDA WorkGraph"
copyright = "2023, Xing Wang"
author = "Xing Wang"

# The full version, including alpha/beta/rc tags
release = "0.0.1"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "nbsphinx",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# Function to copy HTML files
def copy_html_files(src, dst):
    """
    Copy all .html files from src to dst, maintaining the directory structure.
    """
    src_path = Path(src)
    for html_file in src_path.rglob("*.html"):
        relative_path = html_file.relative_to(src_path)
        destination_file = Path(dst) / relative_path
        destination_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(html_file, destination_file)


def setup(app):
    app.connect(
        "build-finished", lambda app, exception: copy_html_files("source", app.outdir)
    )
