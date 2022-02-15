# -*- coding: utf-8 -*-
#
# TRIQS documentation build configuration file

import sys


sys.path.insert(0, "@CMAKE_CURRENT_SOURCE_DIR@/sphinxext")
sys.path.insert(0, "@CMAKE_CURRENT_SOURCE_DIR@/sphinxext/numpydoc")

exclude_patterns = ['_templates', 'examples']

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.intersphinx',
              'sphinx.ext.doctest',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosummary',
              'sphinx.ext.githubpages',
              'sphinx_autorun',
              'myst_parser',
              'nbsphinx',
              'IPython.sphinxext.ipython_console_highlighting',
              'matplotlib.sphinxext.plot_directive',
              'numpydoc']

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

autodoc_default_options = {
    'members': None, # Include all members (methods).
    'special-members': None,
    'exclude-members': '__dict__,__weakref__',
    'undoc-members' : None,
    }

source_suffix = '.rst'

autosummary_generate = True  # Turn on sphinx.ext.autosummary

project = '@PROJECT_NAME@'
version = '@PROJECT_VERSION@'

copyright = 'Copyright (C) 2018-2020, ETH Zurich Copyright (C) 2021, The Simons Foundation authors: A. Hampel, M. Merkel, and S. Beck'

mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=default"
templates_path = ['@CMAKE_CURRENT_SOURCE_DIR@/_templates']

html_theme = 'sphinx_rtd_theme'
html_style = 'css/custom.css'

html_favicon = '@CMAKE_CURRENT_SOURCE_DIR@/logos/favicon.ico'
html_logo = '@CMAKE_CURRENT_SOURCE_DIR@/logos/logo_no_text.png'

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#7E588A',
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

html_show_sphinx = True
html_context = {'header_title': '@PROJECT_NAME@'}
html_static_path = ['@CMAKE_CURRENT_SOURCE_DIR@/_static']
html_sidebars = {'index': ['sideb.html', 'searchbox.html']}

htmlhelp_basename = '@PROJECT_NAME@doc'

intersphinx_mapping = {'python': ('https://docs.python.org/3.8', None), 'triqslibs': ('https://triqs.github.io/triqs/latest', None)}


# open links in new tab instead
from sphinx.writers.html import HTMLTranslator
from docutils import nodes
from docutils.nodes import Element

class PatchedHTMLTranslator(HTMLTranslator):

    def visit_reference(self, node: Element) -> None:
        atts = {'class': 'reference'}
        if node.get('internal') or 'refuri' not in node:
            atts['class'] += ' internal'
        else:
            atts['class'] += ' external'
            # ---------------------------------------------------------
            # Customize behavior (open in new tab, secure linking site)
            atts['target'] = '_blank'
            atts['rel'] = 'noopener noreferrer'
            # ---------------------------------------------------------
        if 'refuri' in node:
            atts['href'] = node['refuri'] or '#'
            if self.settings.cloak_email_addresses and atts['href'].startswith('mailto:'):
                atts['href'] = self.cloak_mailto(atts['href'])
                self.in_mailto = True
        else:
            assert 'refid' in node, \
                   'References must have "refuri" or "refid" attribute.'
            atts['href'] = '#' + node['refid']
        if not isinstance(node.parent, nodes.TextElement):
            assert len(node) == 1 and isinstance(node[0], nodes.image)
            atts['class'] += ' image-reference'
        if 'reftitle' in node:
            atts['title'] = node['reftitle']
        if 'target' in node:
            atts['target'] = node['target']
        self.body.append(self.starttag(node, 'a', '', **atts))

        if node.get('secnumber'):
            self.body.append(('%s' + self.secnumber_suffix) %
                             '.'.join(map(str, node['secnumber'])))

def setup(app):
    app.set_translator('html', PatchedHTMLTranslator)
