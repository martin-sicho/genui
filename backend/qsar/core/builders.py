"""
builders

Created by: Martin Sicho
On: 15-01-20, 12:55
"""
from django.core.exceptions import ImproperlyConfigured

import modelling.core.bases
import modelling.models
from qsar.models import ModelActivity
from . import bases
from qsar import models
from sklearn.model_selection import KFold, StratifiedKFold

class BasicQSARModelBuilder(bases.QSARModelBuilder):

    def build(self) -> models.QSARModel:
        if not self.validation:
            raise ImproperlyConfigured("You cannot build a QSAR model with a missing validation strategy.")

        if not self.molsets:
            raise ImproperlyConfigured("You cannot build a QSAR model without an associated molecule set.")

        self.progressStages = [
            "Fetching activities...",
            "Calculating descriptors..."
        ]
        self.progressStages.extend([f"CV fold {x+1}" for x in range(self.validation.cvFolds)])
        self.progressStages.extend(["Fitting model on the training set...", "Validating on test set..."])
        self.progressStages.extend(["Fitting the final model..."])

        self.recordProgress()
        mols = self.saveActivities()[1]

        self.recordProgress()
        self.calculateDescriptors(mols)

        X_valid = self.X.sample(frac=self.validation.validSetSize)
        X_train = self.X.drop(X_valid.index)
        y_valid = self.y[X_valid.index]
        y_train = self.y.drop(y_valid.index)

        is_regression = self.training.mode.name == modelling.core.bases.Algorithm.REGRESSION
        if is_regression:
            folds = KFold(self.validation.cvFolds).split(X_train)
        else:
            folds = StratifiedKFold(self.validation.cvFolds).split(X_train, y_train)
        for i, (trained, validated) in enumerate(folds):
            self.recordProgress()
            self.fitAndValidate(X_train.iloc[trained], y_train.iloc[trained], X_train.iloc[validated], y_train.iloc[validated], perfClass=modelling.models.ModelPerformanceCV, fold=i + 1)

        model = self.algorithmClass(self)
        self.recordProgress()
        model.fit(X_train, y_train)
        self.recordProgress()
        self.fitAndValidate(X_train, y_train, X_valid, y_valid)
        self.recordProgress()
        return super().build()

    def populateActivitySet(self, aset : models.ModelActivitySet):
        aset.activities.all().delete()
        molecules = aset.molecules.molecules.all()
        self.calculateDescriptors(molecules)
        predictions = self.model.predict(self.getX())

        for mol, prediction in zip(molecules, predictions):
            ModelActivity.objects.create(
                value=prediction,
                type=self.training.modelledActivityType,
                units=self.training.modelledActivityUnits,
                source=aset,
                molecule=mol,
            )

        return aset.activities.all()

