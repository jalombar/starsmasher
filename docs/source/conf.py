# Configuration for the StarSmasher documentation build.
project   = 'StarSmasher'
copyright = '2026, James C. Lombardi Jr. and contributors'
author    = 'James C. Lombardi Jr. and contributors'
release   = 'master'

extensions = [
    'myst_parser',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
    'sphinx_design',
]

myst_enable_extensions = ['deflist', 'colon_fence', 'substitution']

templates_path   = ['_templates']
exclude_patterns = []
source_suffix    = {'.rst': 'restructuredtext', '.md': 'markdown'}

html_theme       = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files   = ['custom.css']
html_title       = 'StarSmasher'
# The logo is optional so the docs still build before the artwork is added.
import os as _os
_logo = _os.path.join(_os.path.dirname(__file__), '_static', 'starsmasher-logo.png')
html_logo = '_static/starsmasher-logo.png' if _os.path.exists(_logo) else None
html_favicon     = None

html_theme_options = {
    'logo_only': True,
    'navigation_depth': 3,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'style_external_links': True,
}

html_context = {
    'display_github': True,
    'github_user': 'jalombar',
    'github_repo': 'starsmasher',
    'github_version': 'master',
    'conf_py_path': '/docs/source/',
}
