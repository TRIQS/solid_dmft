{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}
    :members:

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

{% block classes %}
{% if classes %}
.. rubric:: {{ _('Classes') }}

.. autosummary::
    :toctree:                                         
    :template: autosummary_class_template.rst
    {% for item in classes %}
    {{ item }}
    {%- endfor %}
{% endif %}
{% endblock %}

{% block exceptions %}
{% if exceptions %}
.. rubric:: {{ _('Exceptions') }}

.. autosummary::
    :toctree:                                       
    {% for item in exceptions %}
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

