..  _dev-guide-generators-top:

Molecular Generators
====================

In this tutorial, we will finally focus on the central part
of the GenUI platform. However, it is recommended that you still check the
:ref:`dev-guide-create-qsar-ext` and :ref:`dev-guide-create-compounds-ext`
tutorials first because we will expand upon some of the concepts
presented there.

A lot of the most promising contemporary generators are based on generative
deep neural networks, which are a modern class of machine learning
algorithms. Therefore, in this tutorial we will mainly be concerned
with the integration of such an approach. However, it should be noted that
any generator can be added to the platform even if it is not strictly
based on machine learning. You can skip to chapter :ref:`dev-guide-generators-useit`
if your generator is not based on a machine learning algorithm or does not require
training. However, you might still find :ref:`dev-guide-generators-urls`, :ref:`dev-guide-generators-views` and :ref:`dev-guide-generators-models-serializers` useful if you
are not familiar with Django and Django REST frameworks.

We will not describe the development of a generative extension in exhausting detail here, but you should regard this tutorial more as a case study which highlights important GenUI concepts and features. The subject of our case study will be the
`genui.generators.extensions.genuidrugex` package, which integrates a deep learning
generative algorithm called `DrugEx <https://github.com/XuhanLiu/DrugEx>`_. This approach is based on training a recurrent neural network with a reinforcement learning loop
where a QSAR model provides the environment for the agent. It is implemented on
top of the PyTorch machine
learning library which can take advantage of GPU hardware. Therefore, DrugEx is
a good non-trivial example of a contemporary molecular generator that might be
of interest to users seeking to use the GenUI platform to discover new chemistry.

Just like any other extension, `genuidrugex`
is also a Django application so it has all the expected modules
and functionality that you are probably already familiar with if
you read the previous tutorials. However,
in order to achieve our goals, we will have to do a little bit more work
and customization than we saw previously. Most of this work will focus on creating a customized *model builder* (see :ref:`dev-guide-generators-builder`), which
has a similar purpose as the instances of `MolSetInitializer` that we
saw in :ref:`dev-guide-create-compounds-initializers`. However, we will also
learn more about other features of GenUI that we have not described, yet.

..  _dev-guide-generators-urls:

Defining URLs
-------------

We will take a look at the `genuidrugex.urls`
module first since it clearly showcases the components of the extension
that we will be dealing with in our case study. The file is not long so we can show it here:

..  literalinclude:: ../../../../src/genui/generators/extensions/genuidrugex/urls.py

Before training the generator itself (the *agent*), the `DrugEx approach <https://doi.org/10.1186/s13321-019-0355-6>`_ requires
that an *exploration* and *exploitation* networks are build first. That is why we see
two viewsets (`DrugExNetViewSet` and `DrugExAgentViewSet`) registered for the router in the above code. Each viewset handles
creation and training of either one of the networks or the agent itself.

The :code:`routes` list is simply a list of URLs for API endpoints
that we would like to have available under the root URLs :code:`/generators/drugex/networks/` and :code:`/generators/drugex/agents/`.
You can see that we are using some class views that come from the `genui.models`
application:

    1. `ModelTasksView`: Shows all or only started and running asynchronous Celery tasks attached to an instance of :class:`~genui.models.models.Model` (determined by the :code:`pk` URL argument).
    2. `ModelPerformanceListView`: Shows machine learning model performance for a specified instance of :class:`~genui.models.models.Model`. See :ref:`dev-guide-qsar-metrics` for more info.
    3. `ModelFileView`: :class:`~genui.models.models.Model` instances can have files attached to it. This endpoint handles them.

You do not have to have these views if you do not need them, but it is good to know
about them so that you do not have to implement your own.

..  _dev-guide-generators-views:

Defining Views
--------------

Lets take a closer look at `DrugExNetViewSet` and `DrugExAgentViewSet` in `genuidrugex.views`:

..  literalinclude:: ../../../../src/genui/generators/extensions/genuidrugex/views.py

GenUI already has a viewset defined for model creation, `ModelViewSet`. Actually, the same viewset is used for QSAR modelling
as well (see `genui.qsar.views`) and it turns out
that we do not actually need to modify it much. There is just a few class attributes
that we have to set:

    1.  :code:`queryset`: Defines the objects shown in the list and their order.
    2.  :code:`serializer_class` and :code:`init_serializer_class`: These are serializer
        classes used to represent and create new instances. They can be identical, but
        in the case of models they are often slightly different. Hence, the distinction
        and two parameters.
    3.  :code:`builder_class`: This is a specific attribute for models and is similar
        to the :code:`initializer_class` found in `genui.compounds` (see
        :ref:`dev-guide-create-compounds-initializers`). We will tak about in more
        detail `later <dev-guide-generators-builder>`_.
    4.  :code:`build_task`: We have a specific Celery task defined for DrugEx (:py:func:`~genui.generators.extensions.genuidrugex.tasks.buildDrugExModel()` in `genuidrugex.tasks`).
        Mainly because this extension depends on GPU hardware and as such the build tasks should be submitted
        to the :code:`gpu` queue, which should be consumed by workers with GPUs.

..  _dev-guide-generators-models-serializers:

Models & Serializers
--------------------

Since we are implementing a completely new kind of model,
we have to define the models and serializers we will need.
This part is very specific to the DrugEx extension,
but there are a few model and serializer classes
in the `genui.models` application that you should be
aware of when attempting something similar.

Model Classes
~~~~~~~~~~~~~

:py:class:`~genui.models.models.Model`
______________________________________

A polymorphic model class used to save data
about a machine learning model. In `genuidrugex.models`, we
can add data that is saved to this class by creating a subclass.
For example, the final DrugEx agent can be described like so:

..  code-block:: python

    class DrugExAgent(Model):
        environment = models.ForeignKey(QSARModel, on_delete=models.CASCADE, null=False, related_name='drugexEnviron')
        explorationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExplore')
        exploitationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExploit')

        def getGenerator(self):
            return self.generator.all().get() if self.generator.all().exists() else None

Here we add foreign keys to the `explorationNet` and `exploitationNet` models (represented by `DrugExNet`, also a subclass of :class:`~genui.models.models.Model`). We also
add the `environment`, which is used for the policy gradient calculation in the reinforcement learning loop of the DrugEx algorithm. It is simply a given QSAR model, also defined
as a subclass of :class:`~genui.models.models.Model` (`QSARModel`).

:class:`~genui.models.models.TrainingStrategy` and :class:`~genui.models.models.ValidationStrategy`
___________________________________________________________________________________________________

When a model is created in GenUI several parameters have to be specified
depending on whether the model is to be trained and validated by GenUI or imported
from an external source. Two most important :py:class:`~genui.models.models.Model`
attributes are :py:attr:`~genui.models.models.Model.trainingStrategy` and
:py:attr:`~genui.models.models.Model.validationStrategy`. These should point to database
model instances of :class:`~genui.models.models.TrainingStrategy` and
:class:`~genui.models.models.ValidationStrategy`.

:class:`~genui.models.models.TrainingStrategy` contains data
required for training of a model. Most importantly this is information about the algorithm used and the chosen training parameters, which are tied to this instance with the `ModelParameterValue` model.

On the other hand, :class:`~genui.models.models.ValidationStrategy`
defines how the model should be validated after it is trained. It holds information
about the performance metrics that should be calculated
during validation (see `ModelPerformanceMetric`).

Just like the :class:`~genui.models.models.Model` class, :class:`~genui.models.models.TrainingStrategy` and
:class:`~genui.models.models.ValidationStrategy` are polymorphic as well and you can add new information by subclassing them.

For example, `BasicValidationStrategy` saves information about the
size of the external validation set chosen randomly from the training
data (defined as a fraction of instances) and the number of folds in cross-validation. It is defined
in `genui.models.models` as follows:

..  code-block:: python

    class CV(ValidationStrategy):
        cvFolds = models.IntegerField(blank=False)

        class Meta:
            abstract = True


    class ValidationSet(ValidationStrategy):
        validSetSize = models.FloatField(blank=False) # as

        class Meta:
            abstract = True


    class BasicValidationStrategy(ValidationSet, CV):
        pass

In the `genuidrugex` extension we have the following custom validation
strategy:

..  code-block:: python

    class DrugExValidationStrategy(ValidationStrategy):
        validSetSize = models.IntegerField(default=512, null=True)

It simply defines the size of the test set that will be used to measure
performance after processing one batch of data.

Serializer Classes
~~~~~~~~~~~~~~~~~~

Classes defined in `genuidrugex.serializers` just describe how the models in `genuidrugex.serializers` are transformed to JSON format. This is
a little bit more involved and you are encouraged to study the source code of this
module more closely. We will just follow up on the example above and show how the serializer looks like for `DrugExValidationStrategy`:

..  code-block:: python

    class DrugExValidationStrategySerializer(ValidationStrategySerializer):
        """
        This is used for GET requests to convert a DrugExValidationStrategy model
        to JSON.
        """


        class Meta:
            model = models.DrugExValidationStrategy
            fields = ValidationStrategySerializer.Meta.fields + ("validSetSize",)

    class DrugExValidationStrategyInitSerializer(DrugExValidationStrategySerializer):
        """
        This is the serializer used for POST requests when creating a new instance of
        DrugExValidationStrategy.
        """

        metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=ModelPerformanceMetric.objects.all(), required=False)

        class Meta:
            model = models.DrugExValidationStrategy
            fields = DrugExValidationStrategySerializer.Meta.fields + ("validSetSize",)

Many serializers are already defined in `genui.models.serializers` for
`Model`, `TrainingStrategy`, `ValidationStrategy` (`ModelSerializer`, `TrainingStrategySerializer`, `ValidationStrategySerializer`) and other classes related to machine learning models. You should always derive from these when building
a customized model.

..  _dev-guide-generators-builder:

Implementing a Model Builder
----------------------------

Now it is time to cover the most essential attribute of `ModelViewSet`, :code:`builder_class`. Each GenUI model needs a builder, which is initialized when the Celery build task is executed. QSAR models have their own builder (`QSARModelBuilder`) and chemical space
maps are also models with their own `MapBuilder`. In this chapter, we will cover
how builder can be implemented for a molecular generator.

Lets just briefly remind ourselves of the views we defined earlier:

..  literalinclude:: ../../../../src/genui/generators/extensions/genuidrugex/views.py

In this code, we should now understand everything except the :code:`builder_class` and :code:`build_task`
attributes and the :meth:`~genui.models.views.ModelViewSet.get_builder_kwargs` method.

When the create view of the `DrugExNetViewSet` or `DrugExAgentViewSet` receives a POST request with
data defining a new DrugEx model, the POST data is validated according to the
appropriate serializers and the :code:`create` method of `DrugExNetInitSerializer` or
`DrugExAgentInitSerializer` is called. This method is responsible for the creation of
a new :class:`~genui.generators.extensions.genuidrugex.models.DrugExNet` or
:class:`~genui.generators.extensions.genuidrugex.models.DrugExAgent` instance from the
supplied data. An asynchronous Celery task (:func:`~genui.generators.extensions.genuidrugex.tasks.buildDrugExModel`
defined in `genuidrugex.tasks`) is then added to queue and passed the ID of the created instance. The fully qualified name of the
builder class and the model class (as specified by the :meth:`~genui.models.views.ModelViewSet.get_builder_kwargs` method of `ModelViewSet`)
are passed as well.

We can look at the code of :func:`~genui.generators.extensions.genuidrugex.tasks.buildDrugExModel` to explain what happens once the task is eligible for execution
by a Celery worker:

..  literalinclude:: ../../../../src/genui/generators/extensions/genuidrugex/tasks.py

Once the task is eligible for execution, the function above is executed on the worker and the
information outlined above is passed to it as parameters. You should always make
sure that whatever information you pass from :meth:`~genui.models.views.ModelViewSet.get_builder_kwargs` is serializable as JSON because that is how this information
is passed from the server to the worker.

..  note:: We also specify :code:`gpu` as the queue for this task. This indicates
    that this task prefers to be executed on a worker node with access to GPUs.

Therefore, looking at the code above this is roughly what happens on the worker:

    1. We import the model class from `genuidrugex.models`.
    2. We use :func:`getObjectAndModuleFromFullName` from `genui.utils.inspection`
       to import the correct builder.
    3. We create an instance of the correct builder (`DrugExNetBuilder` or `DrugExAgentBuilder`).
    4. We run the :command:`build` method of the builder to build the model.

Therefore, in order to build models we have to implement a builder, which means
subclassing the :py:class:`~genui.models.genuimodels.bases.ModelBuilder` abstract class and adding the necessary methods. You can find the DrugEx model builders in `genuidrugex.genuimodels.builders`, which is a standard location recognized by :code:`genuisetup`.
You should always define your builders in :code:`{your_extension}.genuimodels.builders`.

Note that there are a few mix-in classes defined in `genui.models.genuimodels.bases` that you can take advantage of when creating your own model builders. For example, in the DrugEx extension we use `ProgressMixIn` that adds methods and attributes to record task progress stages more easily.

The purpose of a model builder is to
prepare model training data with the :meth:`~genui.models.genuimodels.bases.ModelBuilder.getX` and :meth:`~genui.models.genuimodels.bases.ModelBuilder.getY` methods. Values returned by these methods are then fed to the implementation of :meth:`~genui.models.genuimodels.bases.Algorithm.fit` of an :class:`~genui.models.genuimodels.bases.Algorithm` instance (accessible from the :code:`model` attribute of :py:class:`~genui.models.genuimodels.bases.ModelBuilder`). This is exactly what the
reference implementation of :meth:`~genui.models.genuimodels.bases.ModelBuilder.build`
in :class:`~genui.models.genuimodels.bases.ModelBuilder` does:

..  code-block:: python

    def build(self) -> models.Model:
        """
        Build method of the ModelBuilder abstract class.
        """

        self.model.fit(self.getX(), self.getY())
        self.saveFile() # calls self.model.serialize(path_to_model_snapshot)
        return self.instance

In the case of `BasicQSARModelBuilder`, :py:meth:`~genui.qsar.genuimodels.builders.BasicQSARModelBuilder.getX` returns the matrix of
descriptor values for each compound and :py:meth:`~genui.qsar.genuimodels.builders.BasicQSARModelBuilder.getY` returns the class labels. In the DrugEx extension,
we return the parsed data set corpus from :py:meth:`~genui.generators.extensions.genuidrugex.genuimodels.builders.DrugExNetBuilder.getX`.

Also note the `saveFile` method of the builder. This ensures that the fitted
model is serialized on disk after training. Each subclass of :class:`~genui.models.genuimodels.bases.Algorithm` is responsible for implementing the
correct :meth:`~genui.models.genuimodels.bases.Algorithm.serialize` method
so that a proper save file can be generated.

You can see the defined DrugEx algorithms by exploring
the `genuidrugex.genuimodels.algorithms` module. You will notice that
the algorithms implement a few more methods in addition to those described
in our discussion of QSAR model algorithms (see :ref:`dev-guide-create-qsar-ext`).
We have to define new file formats for model serialization with the `DrugExAlgorithm.getFileFormats()` method,
make sure that our algorithm lists only as a generator by providing the correct
mode with `DrugExAlgorithm.getModes()` and we also have to provide a new serializer and deserializer (`DrugExAlgorithm.getSerializer()` and `DrugExAlgorithm.getDeserializer()`)
for the state of the model (used by `Algorithm.serialize()`).

..  _dev-guide-generators-useit:

Using the Trained Generator
---------------------------

So far we have described a possible implementation of a machine learning based generator on the example of the DrugEx extension,
but we have not yet described how to use the trained generator for the creation of new compound sets. You might have
noticed the definition of the `DrugEx` class in `genuidrugex.models`:

..  code-block:: python

    class DrugEx(Generator):
        agent = models.ForeignKey(Model, on_delete=models.CASCADE, null=False, related_name="generator")

        def get(self, n_samples):
            import genui.generators.extensions.genuidrugex.genuimodels.builders as builders
            builder_class = getattr(builders, self.agent.builder.name)
            builder = builder_class(Model.objects.get(pk=self.agent.id))
            samples, valids = builder.sample(n_samples)
            return [x for idx, x in enumerate(samples) if bool(valids[idx])]

We have not discussed this Django model yet because it is not central to our discussion
of training and saving the DrugEx networks, but an instance of this class is
created whenever an instance of `DrugExNet` or :py:class:`~genui.generators.extensions.genuidrugex.models.DrugExAgent` is saved
and in fact it is all we need to register a new generator with GenUI and
generate new compounds with it.

Looking at the implementation above, we see that
the definition of a generator in GenUI is
quite general and straightforward. It is simply any instance of the `Generator`
class. It should always implement the :py:meth:`~genui.generators.models.Generator.get` method which takes only one
argument, the maximum number of compounds to generate. In the case of DrugEx,
the implementation of :py:meth:`~genui.generators.extensions.genuidrugex.models.DrugEx.get` we see above is only a question of importing the correct
DrugEx builder, which we equipped with the :py:meth:`~genui.generators.extensions.genuidrugex.genuimodels.algorithms.DrugExAlgorithm.sample` method that can use the
trained neural network to generate compounds.

After implementing your `Generator` Django model, you should be able to see
it as an option when creating new compound sets with the
API endpoints of the `genui.compounds.extensions.generated` extension.
Note that this is really all you need to define
a generator and if you do not require to train a machine learning model,
your Django application could be really simple and reduced to just implementing
this Django model, creating an appropriate serializer and hooking it up with
a simple Django REST Framework view or viewset.

Conclusion
----------

This was a short tour through how the `genuidrugex` extension was created
and hopefully it is a bit more clear how you can integrate your own generator
in GenUI. The `genui.models` package has much more interesting features
than covered here so you are encouraged to check its own documentation
pages for some guidance on how to implement the functionality that you want.
