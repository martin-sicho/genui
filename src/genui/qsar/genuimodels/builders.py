"""
builders

Created by: Martin Sicho
On: 15-01-20, 12:55
"""
# CHANGE: Added import for traceback to handle exceptions
import traceback

import numpy as np
from django.core.exceptions import ImproperlyConfigured
from rdkit import Chem

from genui.models.genuimodels.bases import Algorithm
from genui.models.models import ModelPerformanceCV, ModelFile
from . import bases
# CHANGE: Updated imports to use more specific model references
from genui.models import models as core_models
from genui.qsar import models as qsar_models
from sklearn.model_selection import KFold, StratifiedKFold, train_test_split
from django.core.files.base import ContentFile

class BasicQSARModelBuilder(bases.QSARModelBuilder):
    # CHANGE: Updated __init__ method to handle multiple validation strategies
    def __init__(self, instance: qsar_models.Model, progress=None, onFitCall=None, validations=None):
        super().__init__(instance, progress, onFitCall)
        self.validations = validations if validations and len(validations) > 0 else self.instance.trainingStrategy.validationStrategies.all()

    def build(self) -> qsar_models.QSARModel:
        if not self.validations:
            raise ImproperlyConfigured("You cannot build a QSAR model without validation strategies.")
        # CHANGE: Now checking for the presence of validation strategies instead of a single validation strategy.
        # This ensures that at least one validation strategy is present before building the model.
        if not self.molsets:
            raise ImproperlyConfigured("You cannot build a QSAR model without an associated molecule set.")

        self.progressStages = [
            "Fetching activities...",
            "Calculating descriptors..."
        ]
        # CHANGE: Generate progress stages for each validation strategy
        for validation in self.validations:
            self.progressStages.extend([f"CV fold {x+1}" for x in range(validation.cvFolds)])
        self.progressStages.extend(["Fitting model on the training set...", "Validating on test set..."])
        self.progressStages.extend(["Fitting the final model...", "Saving the model..."])

        self.recordProgress()
        mols = self.saveActivities()[1]

        self.recordProgress()
        self.calculateDescriptors(mols)

        X = self.getX()
        y = self.getY()

        for validation in self.validations:
            if hasattr(validation, 'validSetSize') and hasattr(validation, 'cvFolds'):
                # Calculate the size of the validation set
                valid_set_size = validation.validSetSize if isinstance(validation.validSetSize, float) else 0.2

                # Use train_test_split for each validation strategy
                X_train, X_valid, y_train, y_valid = train_test_split(
                    X, y, 
                    test_size=valid_set_size, 
                    random_state=42  # You can make this configurable if needed
                )

                # Perform cross-validation
                is_regression = self.training.mode.name == Algorithm.REGRESSION
                if is_regression:
                    folds = KFold(validation.cvFolds).split(X_train)
                else:
                    folds = StratifiedKFold(validation.cvFolds).split(X_train, y_train)

                for i, (train_index, test_index) in enumerate(folds):
                    self.recordProgress()
                    self.fitAndValidate(
                        X_train.iloc[train_index], y_train.iloc[train_index],
                        X_train.iloc[test_index], y_train.iloc[test_index],
                        perfClass=core_models.ModelPerformanceCV, 
                        fold=i + 1
                    )

                # Validate on the held-out validation set
                self.fitAndValidate(X_train, y_train, X_valid, y_valid)

        # Final model fitting on all data
        final_model = self.algorithmClass(self)
        final_model.fit(X, y)
        self._model = final_model

        # Final validation (optional, as it's not truly a validation)
        y_predicted = final_model.predict(X)
        for validation in self.validations:
            self.validate(validation, y, y_predicted)

        self.recordProgress()
        self.saveFile()
        return self.instance

    def saveModel(self, model):
        """
        Save the trained model to a file.

        This method handles the creation and updating of the model file associated
        with the current instance. If no model file exists, it creates one. Then,
        it serializes the model to the file.

        Args:
            model: The trained model object to be saved.

        Note:
            - Uses the first file format associated with the training algorithm.
            - Creates a new ModelFile if one doesn't exist for this instance.
            - Serializes the model using the model's own serialization method.
        """
        # Get the first file format associated with the training algorithm
        model_format = self.training.algorithm.fileFormats.all()[0]

        # Create a new ModelFile if one doesn't exist for this instance
        if not self.instance.modelFile:
            ModelFile.create(
                self.instance,
                f'main.{model_format.fileExtension}',
                ContentFile('placeholder'),  # Initial content, will be overwritten
                kind=ModelFile.MAIN,
                note=f'{self.training.algorithm.name}_main'
            )

        # Get the path of the model file
        path = self.instance.modelFile.path

        # Serialize the model to the file
        model.serialize(path)

        # Save the instance to ensure any changes are persisted
        self.instance.save()

    def predictMolecules(self, mols):
        smiles = []
        failed_indices = []
        for idx, mol in enumerate(mols):
            if mol:
                smiles.append(Chem.MolToSmiles(mol) if type(mol) == Chem.Mol else mol)
            else:
                failed_indices.append(idx)

        predictions = [-1] * len(mols)
        if len(failed_indices) == len(mols):
            return np.array(predictions)

        self.calculateDescriptors(smiles)

        real_predictions = list(self.predict(self.getX()))
        for idx,prediction in enumerate(predictions):
            if idx not in failed_indices:
                predictions[idx] = real_predictions.pop(0)
        assert len(real_predictions) == 0
        return np.array(predictions)
    
    # CHANGE: Updated to use qsar_models instead of models
    def populateActivitySet(self, aset : qsar_models.ModelActivitySet):
        if not self.instance.predictionsType:
            raise Exception("The activity type for QSAR model predictions is not specified.")

        aset.activities.all().delete()
        molecules = aset.molecules.molecules.all()
        predictions = self.predictMolecules(molecules)

        for mol, prediction in zip(molecules, predictions):
            qsar_models.ModelActivity.objects.create(
                value=prediction,
                type=self.instance.predictionsType,
                units=self.instance.predictionsUnits,
                source=aset,
                molecule=mol,
            )

        return aset.activities.all()
    
    # CHANGE: Added new method to handle fitting and validation
    def fitAndValidate(self, X_train, y_train, X_valid, y_valid, y_predicted=None, perfClass=core_models.ModelPerformance, *args, **kwargs):
        if not y_predicted:
            model = self.algorithmClass(self)
            model.fit(X_train, y_train)
            y_predicted = model.predict(X_valid)
        for validation in self.validations:
            self.validate(validation, y_valid, y_predicted, perfClass, *args, **kwargs)
            
    # CHANGE: Added new method to handle individual validation
    def validate(self, validation_strategy, y_validated, y_predicted, perfClass=core_models.ModelPerformance, *args, **kwargs):
        if not validation_strategy:
            raise ImproperlyConfigured(f"No validation strategy is set for model: {repr(self.instance)}")
        for metric_class in self.metricClasses:
            try:
                metric_class(self).save(y_validated, y_predicted, perfClass, *args, **kwargs)
            except Exception as exp:
                print("Failed to obtain values for metric: ", metric_class.name)
                self.errors.append(exp)
                traceback.print_exc()
                continue
            