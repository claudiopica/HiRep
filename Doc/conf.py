#!/usr/bin/env python3

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# -- General configuration ------------------------------------------------
extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.imgmath",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
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
user_doc = "user_guide"
dev_doc = "development_manual"

# General information about the project.
project = "HiRep"
copyright = "2022"
author = "Rudy Arthur, Vincent Drach, Martin Hansen, Sofie Martins, Claudio Pica, Antonio Rago, Fernando Romero-López"
# author = "HiRep Developers"


# Custom CSS
html_static_path = ["_static"]
html_css_files = [
    "css/custom.css",
]

# -- Options for HTML output ----------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes"]

# -- Options for LaTeX output ---------------------------------------------
latex_elements = {
    "sphinxsetup": "VerbatimBorderColor={rgb}{0.9,0.9,0.9}, VerbatimColor={rgb}{0.98,0.98,0.98},  HeaderFamily=\sffamily",
    "preamble": r"""
    \pdfimageresolution=144
    \usepackage{booktabs}
    \usepackage{titlesec}
    \renewcommand{\familydefault}{\sfdefault}
    \renewcommand{\baselinestretch}{1.15}
    \renewcommand{\arraystretch}{1.15}
    \titleformat{\chapter}{\fontsize{27}{27}\selectfont}{\thechapter}{1em}{}
    \titleformat{\section}{\fontsize{20}{20}\selectfont}{\thesection}{1em}{}
    \titleformat{\subsection}{\fontsize{16}{16}\selectfont}{\thesubsection}{1em}{}
""",
    "maketitle": r"""
\makeatletter
\thispagestyle{empty}
{
\centering
\hrulefill
\vspace{1cm}\\
{\fontsize{35}{35}\selectfont\@title}
\vspace{18cm}\\
{\centering\large\textit{\@author}}
\vspace{3em}
}
\makeatother
""",
}
latex_documents = [
    (
        user_doc,
        "hirep_user_guide.tex",
        "HiRep User Guide",
        author,
        "manual",
    ),
    (
        dev_doc,
        "hirep_development_manual.tex",
        "HiRep Development Manual",
        author,
        "manual",
    ),
]

# -- Options for manual page output ---------------------------------------
man_pages = [
    (master_doc, "sphinx-example", "sphinx-example Documentation", [author], 1)
]

# -- Options for Texinfo output -------------------------------------------
texinfo_documents = [
    (
        user_doc,
        "hirep_user_guide",
        "HiRep User Guide",
        author,
        "User Guide",
        "Library for Numerical Simulations of Fermions in Higher Representation on the Lattice",
        "Miscellaneous",
    ),
    (
        dev_doc,
        "hirep_development_manual",
        "HiRep Development Manual",
        author,
        "Development Manual",
        "Library for Numerical Simulations of Fermions on Higher Representation on the Lattice",
        "Miscellaneous",
    ),
]

# -- Extension configuration -------------------------------------------------
import subprocess

subprocess.call("make clean", shell=True)
subprocess.call("cd doxygen ; doxygen", shell=True)

breathe_projects = {"HiRep": "doxygen/build/xml/"}
breathe_default_project = "HiRep"
