#!/usr/bin/env python3

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# -- General configuration ------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
]
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "substitution",
]
myst_heading_anchors = 3

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = [".rst", ".md"]

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "HiRep"
copyright = "2022"
author = "HiRep Developers"

# Custom CSS
html_static_path = ["_static"]
html_css_files = [
    "css/custom.css",
]

# -- Options for HTML output ----------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes"]

# -- Options for LaTeX output ---------------------------------------------
latex_elements = {"preamble": r"\pdfimageresolution=144"}
latex_documents = [
    (
        master_doc,
        "hirep-documenation.tex",
        "HiRep Documentation",
        "HiRep Developers",
        "manual",
    )
]

# -- Options for manual page output ---------------------------------------
man_pages = [
    (master_doc, "sphinx-example", "sphinx-example Documentation", [author], 1)
]

# -- Options for Texinfo output -------------------------------------------
texinfo_documents = [
    (
        master_doc,
        "hiRep-documenation",
        "HiRep Documentation",
        author,
        "sphinx-example",
        "One line description of project.",
        "Miscellaneous",
    ),
]
