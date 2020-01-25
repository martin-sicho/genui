"""
algorithms

Created by: Martin Sicho
On: 14-01-20, 10:16
"""

import modelling.models
from commons.helpers import findClassInModule
from modelling.core.bases import Algorithm, ModelBuilder
from qsar import models
import pandas as pd
from pandas import DataFrame, Series


class DescriptorCalculator:
    group_name = None

    def __init__(self, builder):
        self.builder = builder

    def __call__(self, smiles):
        pass

    @classmethod
    def getDjangoModel(cls) -> models.DescriptorGroup:
        if not cls.group_name:
            raise Exception('You have to specify a name for the descriptor group in its class "group_name" property')
        return models.DescriptorGroup.objects.get_or_create(name=cls.group_name)[0]


class QSARModelBuilder(ModelBuilder):

    @staticmethod
    def findDescriptorClass(name):
        from . import descriptors
        return findClassInModule(DescriptorCalculator, descriptors, "group_name", name)

    def __init__(self, instance: models.QSARModel, progress=None, onFitCall=None):
        super().__init__(instance, progress, onFitCall)

        self.molset = self.instance.molset
        self.descriptorClasses = [self.findDescriptorClass(x.name) for x in self.training.descriptors.all()]

        self.X = None
        self.y = None

    def getX(self) -> DataFrame:
        return self.X

    def getY(self) -> Series:
        return self.y

    def calculateDescriptors(self, mols):
        smiles = [x.canonicalSMILES for x in mols]
        if not self.getX():
            self.X = DataFrame()
            for desc_class in self.descriptorClasses:
                calculator = desc_class(self)
                temp = calculator(smiles)
                temp.columns = [f"{desc_class.group_name}_{x}" for x in temp.columns]
                self.X = pd.concat([self.X, temp], axis=1)
            return self.X

    def saveActivities(self):
        if not self.getY():
            activity_set = self.instance.molset.activities.get()
            compounds, activities = activity_set.cleanForModelling()
            activities = Series(activities)
            if self.training.mode.name == Algorithm.CLASSIFICATION:
                activity_thrs = self.training.activityThreshold
                activities = activities.apply(lambda x : 1 if x >= activity_thrs else 0)
            self.y = activities
            return self.y, compounds

    def fitValidate(self) -> models.QSARModel:
        ret = self.fit()
        self.validate(self.model, self.X, self.y)
        return ret

    def validate(self, model, X : DataFrame, y_truth : Series, perfClass=modelling.models.ModelPerformance, **kwargs):
        if self.training.mode.name == Algorithm.CLASSIFICATION:
            predictions = model.predict(X)[:,1]
            predictions = [1 if x >= 0.5 else 0 for x in predictions]
        else:
            predictions = model.predict(X)
        for metric_class in self.metricClasses:
            try:
                performance = metric_class(self)(y_truth, predictions)
                perfClass.objects.create(
                    metric=modelling.models.ModelPerformanceMetric.objects.get(name=metric_class.name),
                    value=performance,
                    model=self.instance,
                    **kwargs
                )
            except Exception as exp:
                print("Failed to obtain values for metric: ", metric_class.name)
                self.errors.append(exp)
                continue