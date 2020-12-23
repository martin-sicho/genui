..  _install-guide-from-source:

Running from Source
~~~~~~~~~~~~~~~~~~~

This guide will be most relevant for you if you want to set up a development server
from a certain code revision, make your own changes and debug them. For deployment in production or a
quick test run of the whole GenUI framework, you might want to look into the
GenUI :ref:`Docker images <install-guide-docker>`.

Getting the Code
================

The GenUI backend application is `available on GitHub <https://github.com/martin-sicho/genui>`_.
More stable revisions of the code can be found in branch `master <https://github.com/martin-sicho/genui/tree/master>`_. Unstable and experimental code is pushed to the `dev/master <https://github.com/martin-sicho/genui/tree/master>`_ branch before it is merged with `master <https://github.com/martin-sicho/genui/tree/master>`_. You can also checkout a particular release using tags.

For example, the following will checkout the repository and switch to the master development branch:

..  code-block:: bash

    git clone --branch dev/master git@github.com:martin-sicho/genui.git

Installing Dependencies
=======================

The GenUI framework has quite a few dependencies some of which are available from the Anaconda Cloud
while some are available to install by :code:`pip`. We recommend initializing an Anaconda environment
first. You can do so using the :code:`environment.yml` file available at the root of the repository:

..  code-block:: bash

    # while at the root of the repo
    conda env create -f environment.yml

After you have the environment set up, you can install the rest of the dependencies via :code:`pip`
using the :code:`requirements.txt` file, also located at the root of the repository:

..  code-block:: bash

    conda activate genui # activate the environment if you haven't already
    pip install -r requirements.txt

Launching the Application Server
================================

The GenUI backend code comes with a fully configured Django project out of the box. Three configurations are available:

#. *debug* -- Configuration with debugging enabled. Intended for local deployment in development.
#. *staging* -- Staging configuration is intended to be deployed on a remote server, but with debugging enabled.
#. *production* -- Deployment in production.

After choosing the desired configuration, you can use the :envvar:`DJANGO_SETTINGS_MODULE`
environment variable to specify the correct module and run the :code:`manage.py` script
located under :code:`src/`. Configuration modules reside in the
`genui.settings` package. So provided you have `genui` on your :envvar:`PYTHONPATH`,
you can run the Django webserver with the desired configuration as follows:

..  code-block:: bash

    conda activate genui # activate the environment if you haven't already
    export DJANGO_SETTINGS_MODULE=genui.settings.debug # use the debug configuration
    python manage.py runserver # run the development server

..  note:: GenUI relies on a PostgreSQL database for data storage. Here,
    we assume that the database server is available to the application on :code:`localhost:5432`.
    For other configurations, you should use the :envvar:`POSTGRES_DB`,
    :envvar:`POSTGRES_USER`, :envvar:`POSTGRES_PASSWORD` and :envvar:`POSTGRES_HOST`
    environment variables to tell the application what database to look for.

Deploying Celery Workers
------------------------

..  todo:: write this


Creating a Configuration Module
===============================

You can also create your own configuration module. In order to help developers in this task,
we defined the `genui.settings.base` and `genui.settings.genuibase` settings modules
that new modules can inherit from. Feel free to check out the source code of `genui.settings.debug` for
an example of how these base modules can
be used to configure a new project.