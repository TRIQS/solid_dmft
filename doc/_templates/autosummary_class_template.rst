{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
    :members:                               
    :show-inheritance:                     
    :inherited-members:                   
    :special-members: __init__
    :noindex:

{% block methods %}
{% if methods %}
.. autosummary::
    :toctree:                             
    {% for item in methods %}
      ~{{ name }}.{{ item }}
    {%- endfor %}
{% endif %}
{% endblock %}

{% block attributes %}
{% if attributes %}
.. rubric:: {{ _('Attributes') }}

.. autosummary::
    :toctree:                             
    {% for item in attributes %}
      ~{{ name }}.{{ item }}
    {%- endfor %}
{% endif %}
{% endblock %}
