
Integrating a Simple Deep Learning QSAR Model with PyTorch
==========================================================

This tutorial discusses an approach to integrate a simple DNN-based QSAR modeling approach to the GenUI backend. In this tutorial, we will use the same API as shown in the basic guide (see :ref:`dev-guide-create-qsar-ext`), but we will need to apply more complex preprocessing to the input data, which requires proper serialization of multiple supporting files. The example shown here implements a basic feed-forward neural network described by Esben Jannik Bjerrum for his blog "Cheminformania" (available `here <https://www.cheminformania.com/building-a-simple-qsar-model-using-a-feed-forward-neural-network-in-pytorch/>`_). In this tutorial, we assume that the :code:`qsarextra` package presented previously (see :ref:`dev-guide-create-qsar-ext`) is already configured and that the class described here is placed in the :code:`qsarextra.genuimodels.algorithms` submodule.


Creating the DNN Class
----------------------

We start with importing our dependencies and creating the :code:`DNN` class that will implement the methods required by the :py:class:`~genui.models.genuimodels.bases.Algorithm` base class:

..  code-block:: python

    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import VarianceThreshold

    import torch
    import torch.nn as nn
    from torch.utils.data import TensorDataset, DataLoader

    class DNN(Algorithm):
        """
        An example integration of a simple feed forward deep neural network with PyTorch.

        The general approach here follows the simple tutorial written by Esben Jannik Bjerrum for his blog "Cheminformania" (https://www.cheminformania.com/building-a-simple-qsar-model-using-a-feed-forward-neural-network-in-pytorch/).
        """

        name = 'DNN'
        parameters = {
            "n_epochs" : {
                "type" : ModelParameter.INTEGER,
                "defaultValue" : 200
            },
            "hidden_size" : {
                "type" : ModelParameter.INTEGER,
                "defaultValue" : 1024
            },
            "dropout_rate" : {
                "type" : ModelParameter.FLOAT,
                "defaultValue" : 0.8
            },
            "learning_rate" : {
                "type" : ModelParameter.FLOAT,
                "defaultValue" : 0.001
            },
            "batch_size" : {
                "type" : ModelParameter.INTEGER,
                "defaultValue" : 256
            },
            "feature_selector_threshold" : {
                "type" : ModelParameter.FLOAT,
                "defaultValue" : 0.05
            }
        }

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            self._model = None
            self.db_model = self.builder.instance
            self.input_scaler = StandardScaler()
            self.output_scaler = StandardScaler()
            self.feature_selector = VarianceThreshold(threshold=self.params['feature_selector_threshold'])

        @property
        def model(self):
            return self._model

        @classmethod
        def getModes(cls):
            return [cls.REGRESSION,]

        def fit(self, X: DataFrame, y: Series):
            pass

        def predict(self, X: DataFrame) -> Series:
            pass

We already know this structure from the previous tutorial. The only thing new here is that this time we decided to also use the :py:meth:`~genui.models.genuimodels.bases.Algorithm.getModes` method of the :py:class:`~genui.models.genuimodels.bases.Algorithm` class to specify that this algorithm should only be used for regression tasks. With a few changes we could, of course, implement the class to support both classification and regression, but we will leave this as an exercise to the reader to keep things more simple in our discussion.

In the class above, we also already implemented the constructor and the :py:attr:`~genui.models.genuimodels.bases.Algorithm.model` property. Notice that in the constructor we already initialize the data structures that will be used for scaling of the input and output of the neural network. We also use the :code:`VarianceThreshold` filter to remove features with low variance (i.e. fingerprint bits that are never or rarely set for the compounds in the training data). These preprocessing classes will use the training data to derive their parameters and, thus, will have to be serialized with the main model file for later use (see :ref:`dev-guide-advanced-dnn-qsar-serialization`).


Let us now focus on the implementation of :py:meth:`~genui.models.genuimodels.bases.Algorithm.fit` and :py:meth:`~genui.models.genuimodels.bases.Algorithm.predict` methods. We will start with :py:meth:`~genui.models.genuimodels.bases.Algorithm.fit`:

..  code-block:: python

    def fit(self, X: DataFrame, y: Series):
        # get device
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        # process input
        X = self.input_scaler.fit_transform(X)
        X = self.feature_selector.fit_transform(X)
        X = torch.tensor(X, device=device).float()

        # process output
        y = self.output_scaler.fit_transform(y.to_numpy().reshape((-1,1)))
        y = torch.tensor(y, device=device).float()

        dataset = TensorDataset(X, y)
        loader = DataLoader(
            dataset=dataset,
            batch_size=self.params['batch_size'],
            shuffle=True
        )

        #Defining the hyperparameters
        input_size = X.size()[-1]     # The input size should fit our fingerprint size
        hidden_size = self.params['hidden_size']   # The size of the hidden layer
        dropout_rate = self.params['dropout_rate']    # The dropout rate
        learning_rate = self.params['learning_rate']  # The learning rate for the optimizer
        output_size = 1        # This is just a single task, so this will be one
        self._model = self.Net(input_size, hidden_size, dropout_rate, output_size)
        self.model.cuda()

        # start training
        criterion = nn.MSELoss()
        optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
        self.model.train() #Ensure the network is in "train" mode with dropouts active
        epochs = 200
        for e in range(epochs):
            running_loss = 0
            for fps, labels in loader:
                # Training pass
                optimizer.zero_grad() # Initialize the gradients, which will be recorded during the forward pass

                output = self.model(fps) #Forward pass of the mini-batch
                loss = criterion(output, labels) #Computing the loss
                loss.backward() # calculate the backward pass
                optimizer.step() # Optimize the weights

                running_loss += loss.item()
            else:
                if e%10 == 0:
                    print("Epoch: %3i Training loss: %0.2F"%(e,(running_loss/len(loader))))

There are a few minor changes if you compare it with the code from the `original tutorial <https://www.cheminformania.com/building-a-simple-qsar-model-using-a-feed-forward-neural-network-in-pytorch/#mobile-menu>`_. One thing you might wonder about is the call to :code:`self.Net`. This is the same network as defined in the original tutorial, but in our case it only becomes an attribute of our :code:`DNN` class and, thus, we access it through :code:`self`. You can see the full definition later in this tutorial (see :ref:`dev-guide-advanced-dnn-qsar-complete-class`). We only do this to keep everything packed under one class, but it might as well reside under a separate module or package.

Another small change is also the fact that we do not use the validation set to obtain validation loss in each epoch. We omitted this to keep the example as simple as possible, but you could easily do that in your implementation. Therefore, in this case we only obtain validation data through k-fold cross-validation and the independent test set if it is specified by the user. In this example, we could also gather the training loss for each epoch by creating instances of `ModelPerfomanceNN`, which is the Django model class that can be used to save data from neural network training, but we decided to omit this also for the sake of simplicity.

Also note that we scale the inputs as well. In the original tutorial, the input was limited to only fingerprints, but in GenUI we might encounter real-valued fingerprints as well, which might be challenging for the network if not scaled. Therefore, in GenUI we scale the input to ensure more consistent representation across various types of descriptors.

Now it is time to implement the :py:meth:`~genui.models.genuimodels.bases.Algorithm.predict` method of our model, which will be much simpler:

..  code-block:: python

    def predict(self, X: DataFrame) -> Series:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        self.model.eval() # set to evaluation mode

        # process input
        X = self.input_scaler.transform(X)
        X = self.feature_selector.transform(X)
        X = torch.tensor(X, device=device).float()

        # make and process predictions
        preds = self.model(X)
        preds = self.output_scaler.inverse_transform(preds.cpu().detach().numpy())
        preds = Series(preds[:,0])

        return Series(preds[:,0])

In order to make predictions, we just switch the model to evaluation mode, transform the inputs and make the predictions, which need to be transformed back to the original activity values before they are returned.

Now we have a complete implementation required to successfully train the model. You should write a simple unit test to see if everything is working:

..  code-block:: python

    def testDNN(self):
        model = self.createTestQSARModel(
            mode = AlgorithmMode.objects.get(name="regression"),
            algorithm = Algorithm.objects.get(name="DNN"),
            parameters={
                "n_epochs" : 200,
                "hidden_size" : 1024,
                "dropout_rate" : 0.8,
                "learning_rate" : 0.001,
                "batch_size" : 16
            },
            descriptors = [
                DescriptorGroup.objects.get(name="MORGANFP"),
                DescriptorGroup.objects.get(name="RDKit_2D")
            ],
            metrics=[
                ModelPerformanceMetric.objects.get(name="R2"),
                ModelPerformanceMetric.objects.get(name="MSE"),
            ]
        )

 Just add this method to the test calls we already created in the previous tutorial (see :ref:`dev-guide-create-qsar-ext-tests`). You should be able to see the training output as well as the JSON returned after the model is created.

..  _dev-guide-advanced-dnn-qsar-serialization:

Model Serialization for Prediction
----------------------------------

So our model can now be trained and returned from the GenUI REST API, but we will need to make a few adjustments in order to ensure that the predictions are meaningful. We have to serialize the input and output scalers and the feature selector. The :py:class:`~genui.models.genuimodels.bases.Algorithm` base class already implements basic serialization of the main model file (the object returned from the :py:attr:`~genui.models.genuimodels.bases.Algorithm.model` property). We will now extend this functionality and include our data structures as well:

..  code-block:: python

    def serialize(self, filename):
        super().serialize(filename)
        self.save_aux_file(self.input_scaler, 'inscaler')
        self.save_aux_file(self.output_scaler, 'outscaler')
        self.save_aux_file(self.feature_selector, 'featselect')

    def save_aux_file(self, object, note):
        model_format = ModelFileFormat.objects.get(fileExtension='.joblib.gz')

        ret = ModelFile.create(
            self.db_model,
            f'{note}{model_format.fileExtension}',
            ContentFile('placeholder'),
            kind=ModelFile.AUXILIARY,
            note=f'{self.name}_{note}'
        )
        path = ret.file.path
        joblib.dump(
            object
            , path
        )
        return ret

The :py:meth:`~genui.models.genuimodels.bases.Algorithm.serialize` method is implemented by the :py:class:`~genui.models.genuimodels.bases.Algorithm` so we simply override it in our implementation. In addition, we create the :code:`save_aux_file` that handled the creation of auxiliary model files for us. Auxiliary files are simply any files that are associated with a model instance, but do not represent the model itself. These files are instances of `ModelFile`, which holds information about the file format, the kind of the model file (main or auxiliary) and also a note, which can be used to distinguish between various auxiliary files (we use it here to distinguish between our scalers and the feature selector).

..  note:: Every file should have a file format with the file extension defined, which is useful to tell apart different file types. For example, if we wanted to determine the serialization/deserialization method automatically.

In this example we are using the `ModelFile.create` method as a shortcut to create our file. Notice that we use a placeholder file in place of an already opened file with the serialized model. This is done to first create the file database entry along with the path to the Django media file storage. We can then use this path to dump the actual serialized object.

..  note:: We are using the :code:`joblib` library to save our data structures as gziped files. This is the default method used in GenUI, but if you want to serialize your model in a different way, it is completely OK. You just need to add the appropriate `ModelFileFormat` before you create the first instance.

In order to make predictions, we also need to change the :py:meth:`~genui.models.genuimodels.bases.Algorithm.deserialize` method:

..  code-block:: python

    def deserialize(self, filename):
        super().deserialize(filename)

        self.input_scaler = self.load_aux_file('inscaler')
        self.output_scaler = self.load_aux_file('outscaler')
        self.feature_selector = self.load_aux_file('featselect')

    def load_aux_file(self, note):
        file = self.db_model.files.get(note=f'{self.name}_{note}')
        return joblib.load(file.path)

This is pretty straightforward and probably requires no commentary. We just load the appropriate files with :code:`joblib` based on the note we saved during serialization.

..  _dev-guide-advanced-dnn-qsar-complete-class:

Finished Class
--------------

After you complete all of the steps above, the class should be fully operational. This is what the finished implementation will look like:

..  code-block:: python

    class DNN(Algorithm):
        """
        An example integration of a simple feed forward deep neural network with PyTorch.

        The general approach here follows the simple tutorial written by Esben Jannik Bjerrum for his blog "Cheminformania" (https://www.cheminformania.com/building-a-simple-qsar-model-using-a-feed-forward-neural-network-in-pytorch/).
        """

        name = 'DNN'
        parameters = {
            "n_epochs" : {
                "type" : ModelParameter.INTEGER,
                "defaultValue" : 200
            },
            "hidden_size" : {
                "type" : ModelParameter.INTEGER,
                "defaultValue" : 1024
            },
            "dropout_rate" : {
                "type" : ModelParameter.FLOAT,
                "defaultValue" : 0.8
            },
            "learning_rate" : {
                "type" : ModelParameter.FLOAT,
                "defaultValue" : 0.001
            },
            "batch_size" : {
                "type" : ModelParameter.INTEGER,
                "defaultValue" : 256
            },
            "feature_selector_threshold" : {
                "type" : ModelParameter.FLOAT,
                "defaultValue" : 0.05
            }
        }

        class Net(nn.Module):
            def __init__(self, input_size, hidden_size, dropout_rate, out_size):
                super().__init__()
                # Three layers and a output layer
                self.fc1 = nn.Linear(input_size, hidden_size)  # 1st Full-Connected Layer
                self.fc2 = nn.Linear(hidden_size, hidden_size)
                self.fc3 = nn.Linear(hidden_size, hidden_size)
                self.fc_out = nn.Linear(hidden_size, out_size) # Output layer
                #Layer normalization for faster training
                self.ln1 = nn.LayerNorm(hidden_size)
                self.ln2 = nn.LayerNorm(hidden_size)
                self.ln3 = nn.LayerNorm(hidden_size)
                #LeakyReLU will be used as the activation function
                self.activation = nn.LeakyReLU()
                #Dropout for regularization
                self.dropout = nn.Dropout(dropout_rate)

            def forward(self, x):# Forward pass: stacking each layer together
                # Fully connected =&amp;gt; Layer Norm =&amp;gt; LeakyReLU =&amp;gt; Dropout times 3
                out = self.fc1(x)
                out = self.ln1(out)
                out = self.activation(out)
                out = self.dropout(out)
                out = self.fc2(out)
                out = self.ln2(out)
                out = self.activation(out)
                out = self.dropout(out)
                out = self.fc3(out)
                out = self.ln3(out)
                out = self.activation(out)
                out = self.dropout(out)
                #Final output layer
                out = self.fc_out(out)
                return out

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            self._model = None
            self.db_model = self.builder.instance
            self.input_scaler = StandardScaler()
            self.output_scaler = StandardScaler()
            self.feature_selector = VarianceThreshold(threshold=self.params['feature_selector_threshold'])

        @property
        def model(self):
            return self._model

        @classmethod
        def getModes(cls):
            return [cls.REGRESSION,]

        def fit(self, X: DataFrame, y: Series):
            # get device
            device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

            # process input
            X = self.input_scaler.fit_transform(X)
            X = self.feature_selector.fit_transform(X)
            X = torch.tensor(X, device=device).float()

            # process output
            y = self.output_scaler.fit_transform(y.to_numpy().reshape((-1,1)))
            y = torch.tensor(y, device=device).float()

            dataset = TensorDataset(X, y)
            loader = DataLoader(
                dataset=dataset,
                batch_size=self.params['batch_size'],
                shuffle=True
            )

            #Defining the hyperparameters
            input_size = X.size()[-1]     # The input size should fit our fingerprint size
            hidden_size = self.params['hidden_size']   # The size of the hidden layer
            dropout_rate = self.params['dropout_rate']    # The dropout rate
            learning_rate = self.params['learning_rate']  # The learning rate for the optimizer
            output_size = 1        # This is just a single task, so this will be one
            self._model = self.Net(input_size, hidden_size, dropout_rate, output_size)
            self.model.cuda()

            # start training
            criterion = nn.MSELoss()
            optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
            self.model.train() #Ensure the network is in "train" mode with dropouts active
            epochs = 200
            for e in range(epochs):
                running_loss = 0
                for fps, labels in loader:
                    # Training pass
                    optimizer.zero_grad() # Initialize the gradients, which will be recorded during the forward pass

                    output = self.model(fps) #Forward pass of the mini-batch
                    loss = criterion(output, labels) #Computing the loss
                    loss.backward() # calculate the backward pass
                    optimizer.step() # Optimize the weights

                    running_loss += loss.item()
                else:
                    if e%10 == 0:
                        print("Epoch: %3i Training loss: %0.2F"%(e,(running_loss/len(loader))))

        def predict(self, X: DataFrame) -> Series:
            device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
            self.model.eval() # set to evaluation mode

            # process input
            X = self.input_scaler.transform(X)
            X = self.feature_selector.transform(X)
            X = torch.tensor(X, device=device).float()

            # make and process predictions
            preds = self.model(X)
            preds = self.output_scaler.inverse_transform(preds.cpu().detach().numpy())

            return Series(preds[:,0])

        def serialize(self, filename):
            super().serialize(filename)
            self.save_aux_file(self.input_scaler, 'inscaler')
            self.save_aux_file(self.output_scaler, 'outscaler')
            self.save_aux_file(self.feature_selector, 'featselect')

        def save_aux_file(self, object, note):
            model_format = ModelFileFormat.objects.get(fileExtension='.joblib.gz')

            ret = ModelFile.create(
                self.db_model,
                f'{note}{model_format.fileExtension}',
                ContentFile('placeholder'),
                kind=ModelFile.AUXILIARY,
                note=f'{self.name}_{note}'
            )
            path = ret.file.path
            joblib.dump(
                object
                , path
            )
            return ret

        def deserialize(self, filename):
            super().deserialize(filename)

            self.input_scaler = self.load_aux_file('inscaler')
            self.output_scaler = self.load_aux_file('outscaler')
            self.feature_selector = self.load_aux_file('featselect')

        def load_aux_file(self, note):
            file = self.db_model.files.get(note=f'{self.name}_{note}')
            return joblib.load(file.path)

You can test both training of the model as well as predictions in your test suite:

..  code-block:: python

    def testDNN(self):
        model = self.createTestQSARModel(
            mode = AlgorithmMode.objects.get(name="regression"),
            algorithm = Algorithm.objects.get(name="DNN"),
            parameters={
                "n_epochs" : 200,
                "hidden_size" : 1024,
                "dropout_rate" : 0.8,
                "learning_rate" : 0.001,
                "batch_size" : 16
            },
            descriptors = [
                DescriptorGroup.objects.get(name="MORGANFP"),
                DescriptorGroup.objects.get(name="RDKit_2D")
            ],
            metrics=[
                ModelPerformanceMetric.objects.get(name="R2"),
                ModelPerformanceMetric.objects.get(name="MSE"),
            ]
        )

        # use the model to make predictions (tests deserialization)
        post_data = {
            "name": f"Predictions using {model.name}",
            "molecules": self.molset.id
        }
        predict_url = reverse('model-predictions', args=[model.id])
        response = self.client.post(predict_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

