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
`genui.settings` package.

..  note:: You can also create your own configuration module. In order to help developers in this task,
    we defined the `genui.settings.base` and `genui.settings.genuibase` settings modules
    that new modules can inherit from. Feel free to check out the source code of `genui.settings.debug` for
    an example of how these base modules can
    be used to configure a new project.

So provided you have `genui` on your :envvar:`PYTHONPATH`,
you can run the Django webserver with the desired configuration as follows:

..  code-block:: bash

    conda activate genui # activate the environment if you haven't already
    cd src/
    export DJANGO_SETTINGS_MODULE=genui.settings.debug # use the debug configuration
    python manage.py migrate # initialize the database
    python manage.py genuisetup # setup genui extension modules
    python manage.py runserver # run the development server

If everything went well, the backend application should now be accessible
from localhost at the port displayed in the output. If the server is
running on port 8000, you can verify this by going to `<http://localhost:8000/api/>`_,
which will display the backend REST API documentation.`

If this is your first time launching the server,
you also need to run the :code:`migrate` and :code:`genuisetup` commands.
The :code:`genuisetup` command inspects the currently installed GenUI extensions and runs their
setup methods. You should run this command each time you install or update an extension.

..  note:: GenUI relies on a PostgreSQL database for data storage. The 'genui.settings.debug'
    configuration file assumes that the database server is available to the application on :code:`localhost:5432`.
    For other configurations, you should use the :envvar:`POSTGRES_DB`,
    :envvar:`POSTGRES_USER`, :envvar:`POSTGRES_PASSWORD` and :envvar:`POSTGRES_HOST`
    environment variables to tell the application what database to look for and what
    credentials to use. Check the source files of
    :code:`genui.settings.stage` and :code:`genui.settings.prod` for details.

Deploying Celery Workers
------------------------

We have just launched the application server, but we also need to setup at least one
Celery worker to consume background tasks. GenUI uses Redis as a message queue
and by default (see `genui.settings.base`) it looks for the Redis
server running at :code:`redis://localhost:6379` with no password. If you want to override this,
you can use the :envvar:`REDIS_HOST` and :envvar:`REDIS_PASSWORD` environment
variables to customize the host (if no password is specified, none will be given
to the server).

Now you can launch the celery worker like so:

..  code-block:: bash

    cd src/
    export DJANGO_SETTINGS_MODULE=genui.settings.debug
    celery worker -A genui -Q celery,gpu --loglevel=info --hostname=debug-worker@%h

If you are launching only one worker, it should consume both the default :code:`celery`
queue and the :code:`gpu` queue.

Running Tests
=============

During development, it might be useful to run unit tests for the Django project.
You can run all tests for the :code:`genui` project with the
following :code:`manage.py` command:

..  code-block:: bash

    export DJANGO_SETTINGS_MODULE=genui.settings.test
    python manage.py test