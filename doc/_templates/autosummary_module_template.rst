{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Module Attributes

   .. autosummary::
      :toctree:                             
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}

.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: autosummary_module_template.rst 
   :recursive:

{% for item in modules %}
    {{ item }}

{%- endfor %}

{% endif %}
{% endblock %}

