..  _install-guide-by-import:

Installing GenUI Backend as a Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have an existing Django project and would like to integrate some
or all GenUI features to it, this is the guide for you. You can also use it for
development of extensions that you would like to keep as separate packages
without directly extending the GenUI source code. This guide assumes that
you are already familiar with how Django projects and applications work.

Before you begin, make sure you have the :code:`genui` package installed:

..  code-block:: bash

    pip install genui

Changing Your Settings Module
=============================

If you want to install the GenUI applications in an existing project,
you have to include them in :code:`INSTALLED_APPS`
and also make sure some important settings are set for GenUI.

If you want to set things up exactly your way and only want to obtain
the minimal set of parameters GenUI applications need to run,
you can make do by just importing the `genui.settings.genuibase` module
definitions in your :code:`settings.py`:

..  code-block:: python

    # in settings.py of your project
    from genui.settings.genuibase import *

This defines some useful variables, including :code:`GENUI_SETTINGS`.
This variable has to be defined in your :code:`settings.py` and it
is a dictionary containing some basic settings. It also
has an automatically generated list of GenUI extensions,
which you can conveniently append directly
to your :code:`INSTALLED_APPS`:

..  code-block:: python

    INSTALLED_APPS = [
        'django.contrib.admin',
        'django.contrib.auth',
        'django.contrib.contenttypes',
        'django.contrib.sessions',
        'django.contrib.messages',
        'django.contrib.staticfiles',
    ] + GENUI_SETTINGS['APPS']

You can modify
:code:`GENUI_SETTINGS` as you see fit. You can add or remove applications and extensions.
This can be useful when you want to use a different application to manage accounts, for example.
You should append any extensions that you define outside the `genui` package
to this list so that the :code:`manage.py genuisetup` command can run their setup methods. You can
see `genui.settings.genuibase` source for more details about the settings it defines.

Note that using this approach, you will have to define many settings yourself such as
settings for the database, celery application and others. You could save yourself some work by
importing settings from `genui.settings.base` instead, which also serves as an example settings
module.

Adding GenUI URLs
=================

Including the URLs defined by the GenUI applications and extensions is easy.
You can just include the `genui.urls` module in your :code:`urls.py`:

..  code-block::

    # urls.py of your project

    from django.urls import path, include

    urlpatterns = [
        path('genui/', include('genui.urls'))
    ]