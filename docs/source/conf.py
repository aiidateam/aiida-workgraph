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

# version = ""

# The master toctree document.
master_doc = "index"


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
    "sphinx_gallery.gen_gallery",
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

gallery_src_relative_dir = (
    "../gallery"  # relative path of the gallery src wrt. sphinx src
)
sphinx_src_autogen_dirs = [
    "autogen",
    "concept/autogen",
    #"tutorial/autogen",
    "howto/autogen",
    "built-in/autogen",
]
# we mimik the structure in the sphinx src directory in the gallery src directory

# path of python scripts that should be executed
gallery_src_dirs = [
    os.path.join(gallery_src_relative_dir, autogen_dir)
    for autogen_dir in sphinx_src_autogen_dirs
]
sphinx_gallery_conf = {
    "filename_pattern": "/*",
    "examples_dirs": gallery_src_dirs,  # in sphinx-gallery doc referred as gallery source
    "gallery_dirs": sphinx_src_autogen_dirs,  # path to where to gallery puts generated files
}

exclude_patterns = []
# ignore in the autogenerated ipynb files to surpress warning
exclude_patterns.extend(
    [
        os.path.join(sphinx_src_autogen_dir, "*ipynb")
        for sphinx_src_autogen_dir in sphinx_src_autogen_dirs
    ]
)

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaste
html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = [
    "css/theme.css",
]

html_theme_options = {
    "source_repository": "https://github.com/superstar54/aiida-workgraph/",
    "source_branch": "main",
    "source_directory": "docs/source",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/superstar54/aiida-workgraph",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47
7.59.4.07.55-.17.55-.38
0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01
1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95
0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68
0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0
3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013
8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
    # "announcement": "<em>Important</em> announcement!",
}

# pygments_style = "colorful"
# pygments_dark_style = "monokai"


# Function to copy HTML files
def copy_html_files(app, exception):
    """
    Copy all .html files from source to build directory, maintaining the directory structure.
    """
    copy_print_info = "Copying HTML files to build directory"
    print()
    print(copy_print_info)
    print(len(copy_print_info) * "=")
    if exception is not None:  # Only copy files if the build succeeded
        print(
            "Build failed, but we still try to copy the HTML files to the build directory"
        )
    try:
        src_path = Path(app.builder.srcdir)
        build_path = Path(app.builder.outdir)

        copy_print_info = f"Copying html files from sphinx src directory {src_path}"
        print()
        print(copy_print_info)
        print(len(copy_print_info) * "-")
        for html_file in src_path.rglob("*.html"):
            relative_path = html_file.relative_to(src_path)
            destination_file = build_path / relative_path
            destination_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(html_file, destination_file)
            print(f"Copy {html_file} to {destination_file}")

        gallery_src_path = Path(app.builder.srcdir / Path(gallery_src_relative_dir))

        copy_print_info = (
            f"Copying html files from gallery src directory {gallery_src_path} to build"
        )
        print()
        print(copy_print_info)
        print(len(copy_print_info) * "-")
        for html_file in gallery_src_path.rglob("*.html"):
            relative_path = html_file.relative_to(gallery_src_path)
            destination_file = build_path / relative_path
            destination_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(html_file, destination_file)
            print(f"Copy {html_file} to {destination_file}")
    except Exception as e:
        print(f"Failed to copy HTML files: {e}")


def setup(app):
    app.connect("build-finished", copy_html_files)
