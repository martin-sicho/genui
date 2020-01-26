"""
builders

Created by: Martin Sicho
On: 15-01-20, 12:55
"""
import modelling.core.bases
import modelling.models
from . import bases
from qsar import models
from sklearn.model_selection import KFold, StratifiedKFold

class BasicQSARModelBuilder(bases.QSARModelBuilder):

    def __init__(
            self
            , *args
            , **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.progressStages = [
            "Fetching activities...",
            "Calculating descriptors..."
        ]
        self.progressStages.extend([f"CV fold {x+1}" for x in range(self.validation.cvFolds)])
        self.progressStages.extend(["Fitting model on the training set...", "Validating on test set..."])
        self.progressStages.extend(["Fitting the final model..."])

    def build(self, callback=None) -> models.QSARModel:
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

