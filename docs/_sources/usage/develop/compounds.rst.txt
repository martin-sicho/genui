..  _dev-guide-create-compounds-ext:

Compounds
=========

The GenUI backend server already provides a few useful
extensions that allow creation of compound sets
from various sources (see `genui.compounds.extensions`).
A *compound set* is an important data structure in GenUI.
It is a way of grouping and organizing compounds
for the purpose of QSAR modelling or training
molecular generators in a GenUI project.
Similarly, there are also *activity sets*
which group biological activities of compounds.

In this tutorial, we will be showing the implementation
of a simple extension that will allow us to
upload compounds and their activity data
in JSON via the REST API. We will call it :code:`jsonimport`.

Creating the Extension
----------------------

The process of creating an extension is no different
from the approach we already outlined :ref:`before <dev-guide-create-extension>`:

..  code-block:: bash

    cd src/
    mkdir genui/compounds/extensions/jsonimport
    python manage.py startapp jsonimport compounds/extensions/jsonimport

Also, do not forget to add your package to the :code:`__init__.py`
of :code:`genui.compounds.extensions`:

..  code-block:: python

    """
    __init__.py in src/genui/compounds/extensions/

    """

    __all__ = ('chembl', 'generated', 'sdf', 'csvimports', 'jsonimport')

..  _dev-guide-create-compounds-urls:

Setting URL Prefix
------------------

In order to be able to upload JSON data to the application, the appropriate
REST API endpoint will need to be set up. We will also likely want endpoints
to do other things with the data after upload (querying/updating/deleting/...).

Luckily, the `genui.compounds` package and the `Django Rest Framework <https://www.django-rest-framework.org/>`_
make this job a little easier. You will only need to define a
URL prefix for the endpoints and attach a customized `viewset <https://www.django-rest-framework.org/api-guide/viewsets/>`_. You can do all
this in the :code:`urls.py` module of your extension app
(create this file if it does not exist):

..  code-block:: python

    """
    urls.py in src/genui/compounds/extensions/jsonimport/

    """
    from django.urls import path, include
    from rest_framework import routers
    from . import views

    router = routers.DefaultRouter()
    router.register(r'sets/json', views.JSONSetViewSet, basename='jsonSet')

    urlpatterns = [
        path('', include(router.urls)),
    ]

Here, we are using a *router*, which sets up the endpoints automatically
for us. All we need to do is implement the controllers in
the :code:`JSONSetViewSet` (see :ref:`dev-guide-create-compounds-viewset`).
Here, we also provide a basename (:code:`jsonSet`) for the routes.
This will allow us to easily get the URL of the appropriate endpoint
with the `reverse <https://docs.djangoproject.com/en/3.1/ref/urlresolvers/#reverse>`_ function in Django.

The last step is to append the router URLs to the :code:`urlpatterns`
variable, which is a special variable Django uses to pick
up URL definitions. The patterns you have defined
here are automatically appended under :code:`/compounds/`
by the `genui.compounds` application (see :ref:`dev-guide-create-compounds-genuisetup`). Therefore, we will find our
endpoints under :code:`/compounds/sets/json/` along the endpoints
for other extensions that import SDF or CSV files.

..  _dev-guide-create-compounds-viewset:

Implementing a Compound Set Viewset
-----------------------------------

Lets now look at how we can implement the :code:`JSONSetViewSet`,
which will power our endpoints. We will define it in
the :code:`views.py` module of our extension:

..  code-block:: python

    """
    views.py in src/genui/compounds/extensions/jsonimport/

    Viewsets of the jsonimport package.
    """


    from genui.compounds.extensions.jsonimport.initializer import JSONSetInitializer
    from genui.compounds.extensions.jsonimport.models import JSONMolSet
    from genui.compounds.extensions.jsonimport.serializers import JSONMolSetSerializer
    from genui.compounds.views import BaseMolSetViewSet


    class JSONSetViewSet(BaseMolSetViewSet):
        queryset = JSONMolSet.objects.all() # a Django model queryset defining this compound set
        serializer_class = JSONMolSetSerializer # JSON object serializer
        initializer_class = JSONSetInitializer # JSON compound set initializer

        def get_initializer_additional_arguments(self, validated_data):
            """
            This method can be used to pass extra arguments
            to the compound set initializer implemented
            by *initializer_class*.

            Parameters
            ----------
            validated_data
                validated data according to the *serializer_class*
            Returns
            -------
            parameters : dict
                keyword parameters for the __init__ method of the *initializer_class*
            """

            return {
                # get the molecules from the validated JSON data
                "molecules" : validated_data["molecules"],
            }

As you can see, we can implement a viewset with only a few lines of code.
The `BaseMolSetViewSet` class from GenUI already handles quite a lot
for us. We just need to tell it a few important things:

    1. *Specify the database query*: In order to save and query the uploaded
    compounds, a `Django model <https://docs.djangoproject.com/en/3.1/topics/db/models/>`_ needs to be created that maps a Python class
    to a database table. The :code:`queryset` parameters
    specifies a Django `queryset <https://docs.djangoproject.com/en/3.0/ref/models/querysets/>`_ that will be used to get the
    model instances for this viewset. We will cover this in more detail later: :ref:`dev-guide-create-compounds-models`.

    2. *Define a serializer class*: `Serializers <https://www.django-rest-framework.org/api-guide/serializers/>`_ are objects
    that can map Django models to JSON objects and vice versa. We will need to define a serializer for the :code:`JSONMolSet` Django model in
    :ref:`dev-guide-create-compounds-serializers` and specify it here.

    3. *Define an initializer class*: An initializer is a
    concept coming from the GenUI framework. An object of this class
    handles the creation of a compound set from the uploaded compounds.
    We will show how to implement it in our case later in this tutorial:
    :ref:`dev-guide-create-compounds-initializers`.

..  _dev-guide-create-compounds-models:

Defining Django Models
----------------------

The GenUI framework already has defined data structures
for storage of chemical data. All compounds are
saved as defined by the `Molecule` Django model
class. This model is polymorphic so you can
subclass it and add your own database fields.
You can do it with the `MolSet` model just
as well. For the purpose of this tutorial,
these two classes are all we will need:

..  code-block:: python

    """
    models.py in src/genui/compounds/extensions/jsonimport/
    """

    from django.db import models
    from genui.compounds.models import Molecule, MolSet


    class JSONMolecule(Molecule):
        name = models.CharField(blank=True, null=False, max_length=1024)

    class JSONMolSet(MolSet):
        pass

In this simple case we did not really change the implementation of `MolSet`
in :code:`JSONMolSet`.
However, it is still a good idea to create a separate model since we
might want to extend it in the future and it is also an easy way
to track data uploaded with our extension.

In the case of the `Molecule` model, we only add one field for the name
of the compound in :code:`JSONMolecule`.

..  _dev-guide-create-compounds-serializers:

Defining Serializers
--------------------

Now that we have our models, it is time to tell Django how it should
translate them to the JSON format. That is the purpose of serializers
and GenUI has the `MolSetSerializer`, `MoleculeSerializer` and `ActivitySerializer` classes already implemented. All we need is
to customize them to fit our models:

..  code-block:: python

    """
    serializers.py in src/genui/compounds/extensions/jsonimport

    """
    from rest_framework import serializers

    from genui.compounds.extensions.jsonimport.models import JSONMolecule, JSONMolSet
    from genui.compounds.models import Activity
    from genui.compounds.serializers import MolSetSerializer, MoleculeSerializer, ActivitySerializer

    class JSONActivitySerializer(ActivitySerializer):
        """
        A simplified serializer for activity. We only
        use three fields from the parent:

            1. the numerical value of the activity
            2. the type of the activity (i.e. IC50)
            3. the units of activity for this value (i.e. nM)
        """

        class Meta:
            model = Activity
            fields = ('value', 'type', 'units')

    class JSONMolSerializer(MoleculeSerializer):
        """
        A simplified serializer for compounds.

        Note that we do not need to specify the name field explicitly.
        The framework picks it up automatically from the *JSONMolecule* model.

        We also serialize the activities as a list of *JSONActivitySerializer*
        instances.
        """

        smiles = serializers.CharField(required=True)
        activities = JSONActivitySerializer(many=True)

        class Meta:
            model = JSONMolecule
            fields = ('id', 'name', 'smiles', 'activities')

    class JSONMolSetSerializer(MolSetSerializer):
        """
        A compound set needs to have a few fields specified
        for successful creation. So in this case we take
        them from the *MolSetSerializer* explicitly and
        also add a list of molecules as specified by
        *JSONMolSerializer*.
        """

        molecules = JSONMolSerializer(many=True)

        class Meta:
            model = JSONMolSet
            fields = MolSetSerializer.Meta.fields + ('molecules',)
            read_only_fields = ('created', 'updated')

        def create(self, validated_data):
            """
            Create an instance of JSONMolSet from the validated data.

            Parameters
            ----------
            validated_data : dict
                Validated and parsed data from the JSON object obtained via POST.

            Returns
            -------
                model_instance : JSONMolSet
            """

            ModelClass = self.Meta.model
            return ModelClass.objects.create(
                name=validated_data["name"]
                , description=validated_data["description"]
                , project=validated_data["project"]
            )

This should allow the application to validate and parse the following data, for example:

..  code-block:: python

    {
        "name": "Test JSON Molecule Set",
        "description": "My molecule set for testing...",
        "project": 1, # id of the project to attach this compound set to
        "molecules" : [
            {
                "name": "Vismodegib",
                "smiles": "CS(=O)(=O)C1=CC(=C(C=C1)C(=O)NC2=CC(=C(C=C2)Cl)C3=CC=CC=N3)Cl",
                "activities": []
            },
            {
                "name": "Captopril",
                "smiles": "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O",
                "activities": [
                    {
                        "value": 20.0,
                        "type": {
                            "value": "IC50"
                        },
                        "units": {
                            "value": "nM"
                        }
                    },
                    {
                        "value": 7.7,
                        "type": {
                            "value": "pIC50"
                        },
                        "units": None
                    }
                ]
            },
            {
                "name": "Nimesulide",
                "smiles": "CS(=O)(=O)Nc1ccc([N+](=O)[O-])cc1Oc1ccccc1",
                "activities": [
                    {
                        "value": 11826.0,
                        "type": {
                            "value": "Ki"
                        },
                        "units": {
                            "value": "nM"
                        }
                    },
                    {
                        "value": 4.93,
                        "type": {
                            "value": "pKi"
                        },
                        "units": None
                    }
                ]
            },
        ]
    }

You can see the descriptions and implementations of the `MolSetSerializer`, `MoleculeSerializer` and `ActivitySerializer` classes to get a better
idea of what other fields they define and what purpose they serve. Take a look at the `genui.compounds.serializers` package
for more info about other serializers as well.

..  _dev-guide-create-compounds-initializers:

Defining Compound Set Initializer
---------------------------------

Looking at the :code:`create` method of :code:`JSONMolSetSerializer` above,
we can finally see how a compound set is initialized. However, we do not yet
see how we can add the compounds we have uploaded to it. Populating a compound
set with new compounds is the responsibility of a compound set initializer.

An initializer is any class derived from the `MolSetInitializer` abstract base class. In particular, we need to implement the `MolSetInitializer.populateInstance` and
`MolSetInitializer.updateInstance` methods. In our case, we will also have to change the :code:`__init__` method because we are also passing the :code:`molecules`
from our POST request via the `BaseMolSetViewSet.get_initializer_additional_arguments` of our viewset (see :ref:`dev-guide-create-compounds-viewset`). An example
implementation of the initializer in our simple example could look like this:

..  code-block:: python

    """
    initializer.py in src/genui/compounds/extensions/jsonimport/

    """
    from genui.compounds.extensions.jsonimport.models import JSONMolecule
    from genui.compounds.initializers.base import MolSetInitializer
    from genui.compounds.models import ActivitySet, Activity, ActivityTypes, ActivityUnits


    class JSONSetInitializer(MolSetInitializer):
        """
        Our initializer. It takes a set of molecules
        as it was parsed from the JSON POST request.
        """

        def __init__(self, *args, molecules=tuple(), **kwargs):
            """

            Parameters
            ----------
            args
                positional arguments
            molecules
                as parsed from the JSON request and supplied by *get_initializer_additional_arguments* of the viewset
            kwargs
                any additional keyword arguments we do not care about
            """

            super().__init__(*args, **kwargs) # arguments that we do not want are passed to the base class
            self.molecules = molecules # save the data to be parsed

        def populateInstance(self) -> int:
            """
            Called when a new compound set is created.

            Returns
            -------
            count : int
                number of unique molecules found in the data
            """

            activity_set = None
            for idx, mol_data in enumerate(self.molecules):
                # note current progress
                progress = 100 * idx / len(self.molecules)
                msg = f"Saving molecule: {mol_data['smiles']}"
                self.progress_recorder.set_progress(progress, 100, description=msg)
                print(msg, f"({progress})")

                # create model instance
                mol_instance = self.addMoleculeFromSMILES(
                    mol_data['smiles'],
                    JSONMolecule,
                    {
                        "name" : mol_data['name']
                    }
                )

                # attach activities
                if mol_data['activities']:
                    if not activity_set:
                        activity_set = ActivitySet.objects.create(
                            name=f"{self.instance.name} - imported activities",
                            description="Activities, which were imported with the data.",
                            project=self.instance.project,
                            molecules=self.instance
                        )

                    for activity in mol_data['activities']:
                        type_ = ActivityTypes.objects.get_or_create(value=activity['type']['value'])[0]
                        units = None
                        if activity['units'] and activity['units']['value']:
                            units = ActivityUnits.objects.get_or_create(value=activity['units']['value'])[0]
                        Activity.objects.create(
                            value=activity['value'],
                            units=units,
                            type=type_,
                            source=activity_set,
                            molecule=mol_instance
                        )


            return self.unique_mols

        def updateInstance(self) -> int:
            """
            Called when a compound set is updated.

            For simplicity, we just remove all original data and populate again.
            Returns
            -------
            count : int
                number of compounds changed
            """

            self.instance.activities.all().delete()
            self.instance.molecules.clear()
            return self.populateInstance()

In comparison to the implementations we have seen so far, there is quite
a lot going on, but it actually is no magic. Lets focus on the
`MolSetInitializer.populateInstance` method because it showcases
the most important features.

We begin with looping over the compounds found in the data:

..  code-block:: python

    activity_set = None
    for idx, mol_data in enumerate(self.molecules):
        # note current progress
        progress = 100 * idx / len(self.molecules)
        msg = f"Saving molecule: {mol_data['smiles']}"
        self.progress_recorder.set_progress(progress, 100, description=msg)
        print(msg, f"({progress})")

We use the :code:`progress_recorder` argument to record our progress.
This is important since importing compounds is done asynchronously
inside a Celery task and the progress recorder is used
to propagate task progress data to the GenUI REST API
services reporting on the status and progress of tasks.

Next, we create a :code:`JSONMolecule` instance from the SMILES string
provided in the JSON data:

..  code-block:: python

    # create model instance
    mol_instance = self.addMoleculeFromSMILES(
        mol_data['smiles'],
        JSONMolecule,
        {
            "name" : mol_data['name']
        }
    )

It is important to do so using the `addMoleculeFromSMILES` method.
This method standardizes the structure of the compound using the
`ChEMBL Structure Pipeline <https://github.com/chembl/ChEMBL_Structure_Pipeline>`_ and saves it into the database.
By calling this method you ensure that there is consistency in the way
structures are stored. All you have to do is specify the Django model
class to use and any extra arguments that should be passed to its constructor.

Finally, we attach the activities to the created compounds if any are found
in the data:

..  code-block:: python

    # attach activities
    if mol_data['activities']:
        if not activity_set:
            activity_set = ActivitySet.objects.create(
                name=f"{self.instance.name} - imported activities",
                description="Activities, which were imported with the data.",
                project=self.instance.project,
                molecules=self.instance
            )

        for activity in mol_data['activities']:
            type_ = ActivityTypes.objects.get_or_create(value=activity['type']['value'])[0]
            units = None
            if activity['units'] and activity['units']['value']:
                units = ActivityUnits.objects.get_or_create(value=activity['units']['value'])[0]
            Activity.objects.create(
                value=activity['value'],
                units=units,
                type=type_,
                source=activity_set,
                molecule=mol_instance
            )

Note that we have to create the activity set first. This will ensure
that we can easily distinguish the imported activities from the
activities calculated from a QSAR model, for example. After that
it is all just a matter of creating the `Activity` instances
and supplying the correct data to the :code:`create` method.

..  _dev-guide-create-compounds-genuisetup:

Setting Up the Extension
------------------------

Now we have almost everything in place to put our extension to
use. The only thing left to do is to create :code:`genuisetup.py`:

..  code-block:: python

    """
    genuisetup.py in src/genui/compounds/extensions/jsonimport/

    Created by: Martin Sicho
    On: 1/12/21, 9:49 AM
    """

    PARENT = 'genui.compounds'

    def setup(*args, **kwargs):
        from . import models
        from genui.utils.init import createGroup
        createGroup(
            "GenUI_Users",
            [
                models.JSONMolecule,
                models.JSONMolSet
            ]
        )

The :code:`PARENT` attribute tells GenUI that this extension is meant
to be as a submodule for the `genui.compounds` application and, thus,
all URLs defined in the extension will be prefixed with :code:`/compounds/`.

The `createGroup` function manages user permissions.
Every extension that defines new Django models should specify
permissions for the "GenUI_Users" group. This determines
what API methods will be available to users. Calling
the `createGroup` function like this gives all "GenUI_Users"
read and write permissions for our newly defined model classes.

Finally, we can migrate the database and run the setup :code:`genuisetup` command to make
sure correct permissions are applied:

..  code-block:: bash

    python manage.py makemigrations
    python manage.py migrate
    python manage.py genuisetup

If you are running the server locally, you should now be able to see the appropriate
REST API endpoints documented at :code:`http://localhost:{your_port}/api/`.

Unit Testing
------------

Just like with other types of extensions, it is a good idea to test them. You can use the following
unit test as a template:

..  code-block:: python

    """
    tests.py in src/genui/compounds/extensions/jsonimport

    """
    import json

    from django.urls import reverse
    from rest_framework.test import APITestCase

    from genui.projects.tests import ProjectMixIn


    class ChEMBLMolSetTestCase(ProjectMixIn, APITestCase):
        """
        We use the 'ProjectMixIn' class to automatically get
        a project instance initialized before our tests.
        It will become available from 'self.project'
        """

        def test_json_upload(self):
            post_data = {
                "name": "Test JSON Molecule Set",
                "description": "My molecule set for testing...",
                "project": self.project.id, # get id of our test project
                "molecules" : [
                    {
                        "name": "Vismodegib",
                        "smiles": "CS(=O)(=O)C1=CC(=C(C=C1)C(=O)NC2=CC(=C(C=C2)Cl)C3=CC=CC=N3)Cl",
                        "activities": [] # our serializer should allow empty activities
                    },
                    {
                        "name": "Captopril",
                        "smiles": "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O",
                        "activities": [
                            {
                                "value": 20.0,
                                "type": {
                                    "value": "IC50"
                                },
                                "units": {
                                    "value": "nM"
                                }
                            },
                            {
                                "value": 7.7,
                                "type": {
                                    "value": "pIC50"
                                },
                                "units": None # the serializer should also allow unspecified units
                            }
                        ]
                    },
                    {
                        "name": "Nimesulide",
                        "smiles": "CS(=O)(=O)Nc1ccc([N+](=O)[O-])cc1Oc1ccccc1",
                        "activities": [
                            {
                                "value": 11826.0,
                                "type": {
                                    "value": "Ki"
                                },
                                "units": {
                                    "value": "nM"
                                }
                            },
                            {
                                "value": 4.93,
                                "type": {
                                    "value": "pKi"
                                },
                                "units": None
                            }
                        ]
                    },
                ]
            }

            # create new compound set instance
            response = self.client.post(reverse('jsonSet-list'), post_data, format='json')
            self.assertEqual(response.status_code, 201)
            set_id = response.data['id']
            act_set_id = response.data['activities'][0]

            # use the detail view to fetch the created instance
            response = self.client.get(reverse('jsonSet-detail', args=[set_id]))
            print(json.dumps(response.data, indent=4))
            self.assertEqual(response.status_code, 200)

            # get summary of the uploaded activity data
            summary_url = reverse('activitySet-summary', args=[act_set_id])
            response = self.client.get(summary_url)
            print(json.dumps(response.data, indent=4))
            self.assertEqual(response.status_code, 200)

Note the use of :code:`jsonSet-{identifier}` to denote the proper view in our viewset (see :ref:`dev-guide-create-compounds-urls`). This naming convention is a feature of the `BasicRouter <https://www.django-rest-framework.org/api-guide/routers/#defaultrouter>`_.

