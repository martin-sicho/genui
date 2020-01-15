"""
builders

Created by: Martin Sicho
On: 15-01-20, 12:55
"""
from random import randint

from . import bases
from qsar import models
from sklearn.model_selection import KFold, StratifiedKFold

class BasicQSARModelBuilder(bases.QSARModelBuilder):

    def __init__(
            self
            , training: models.QSARTrainingStrategy
            , validation: models.BasicValidationStrategy
            , onFitCall = None
            , onCVFitCall = None
            , onValidFitCall = None
    ):
        super().__init__(training, validation, onFitCall)
        self.onCVFitCall = onCVFitCall
        self.onValidFitCall = onValidFitCall
        random_state = randint(0, 2**32 - 1)
        self.X_valid = self.X.sample(frac=self.validation.validSetSize, random_state=random_state)
        self.X_train = self.X.drop(self.X_valid.index)
        self.y_valid = self.y.sample(frac=self.validation.validSetSize, random_state=random_state)
        self.y_train = self.y.drop(self.y_valid.index)

    def fitValidate(self) -> models.QSARModel:

        is_regression = self.training.mode == models.TrainingStrategy.REGRESSION
        if is_regression:
            folds = KFold(self.validation.cvFolds).split(self.X_train)
        else:
            folds = StratifiedKFold(self.validation.cvFolds).split(self.X_train, self.y_train)
        for i, (trained, validated) in enumerate(folds):
            model = self.algorithmClass(self.training)
            model.fit(self.X_train[trained], self.y_train[trained],)
            self.validate(model, self.X_train[validated], self.y_train[validated], perfClass=models.ModelPerformanceCV, fold=i)
            self.onCVFitCall(self, i)

        model = self.algorithmClass(self.training)
        model.fit(self.X_train, self.y_train)
        self.validate(model, self.X_valid, self.y_valid)
        self.onValidFitCall(self)
        return self.fit()

