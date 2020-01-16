"""
builders

Created by: Martin Sicho
On: 15-01-20, 12:55
"""
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
        self.stages = [
            "Fetching activities...",
            "Calculating descriptors..."
        ]
        self.stages.extend([f"CV fold {x+1}" for x in range(self.validation.cvFolds)])
        self.stages.extend(["Fitting model on the training set...", "Validating on test set..."])
        self.stages.extend(["Fitting the final model..."])

    def fitValidate(self) -> models.QSARModel:
        self.recordProgress()
        mols = self.saveActivities()[1]

        self.recordProgress()
        self.calculateDescriptors(mols)

        X_valid = self.X.sample(frac=self.validation.validSetSize)
        X_train = self.X.drop(X_valid.index)
        y_valid = self.y[X_valid.index]
        y_train = self.y.drop(y_valid.index)

        is_regression = self.training.mode == bases.Algorithm.REGRESSION
        if is_regression:
            folds = KFold(self.validation.cvFolds).split(X_train)
        else:
            folds = StratifiedKFold(self.validation.cvFolds).split(X_train, y_train)
        for i, (trained, validated) in enumerate(folds):
            self.recordProgress()
            model = self.algorithmClass(self.training)
            model.fit(X_train[trained], y_train[trained],)
            self.validate(model, X_train[validated], y_train[validated], perfClass=models.ModelPerformanceCV, fold=i)

        model = self.algorithmClass(self.training)
        self.recordProgress()
        model.fit(X_train, y_train)
        self.recordProgress()
        self.validate(model, X_valid, y_valid)
        self.recordProgress()
        return self.fit()

