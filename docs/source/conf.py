# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from sphinx_pyproject import SphinxConfig

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
sys.path.insert(0, os.path.abspath("../../src"))
from unfoc.read import __version__ as dynamic_release

# --- Work around sphinx autodoc type-comment signature crash ---
try:
    from sphinx.ext.autodoc import type_comment as _tc

    _real_update = _tc.update_annotations_using_type_comments

    def _safe_update_annotations_using_type_comments(*args, **kwargs):
        try:
            return _real_update(*args, **kwargs)
        except Exception:
            return None

    _tc.update_annotations_using_type_comments = _safe_update_annotations_using_type_comments
except Exception:
    pass
# --------------------------------------------------------------

config = SphinxConfig("../../pyproject.toml", globalns=globals(),
    config_overrides={
        "release": dynamic_release,
        "version": ".".join(dynamic_release.split(".")[:2])  # Optional: e.g., '1.2' from '1.2.3'
    }
)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'numpydoc',
]

add_module_names = False     
autosummary_generate_overwrite = False
autosummary_generate = True
autosummary_context = {"show_docstring": True}

autodoc_default_options = {
    'members': True,
    "undoc-members": False,
    "private-members": False,
    "special-members": "__init__",
    "inherited-members": True,
    "show-inheritance": True,
}

napoleon_numpy_docstring = True
napoleon_google_docstring = False
numpydoc_show_class_members = False
autodoc_typehints = "none"


templates_path = ['_templates']
exclude_patterns = []
html_extra_path = ["../dataflow"]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

