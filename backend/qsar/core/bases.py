"""
core

Created by: Martin Sicho
On: 14-01-20, 10:16
"""

import modelling.models
from commons.helpers import findClassInModule
from compounds.models import Molecule
from modelling.core.bases import Algorithm, CompleteBuilder
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

class DescriptorBuilderMixIn:

    @staticmethod
    def findDescriptorClass(name):
        from . import descriptors
        return findClassInModule(DescriptorCalculator, descriptors, "group_name", name)

    def __init__(self, instance: models.Model, progress=None, onFitCall=None):
        super().__init__(instance, progress, onFitCall)
        self.molset = self.instance.molset
        self.descriptorClasses = [self.findDescriptorClass(x.name) for x in self.training.descriptors.all()]

        self.X = None
        self.y = None

    def calculateDescriptors(self, mols):
        """
        Calculate descriptors for the given molecules
        and save them as X in this instance.

        :param mols: List of molecules to save as X. Can be either instances of Molecule or smiles strings
        :return:
        """

        smiles = [x.canonicalSMILES if isinstance(x, Molecule) else x for x in mols]
        self.X = DataFrame()
        for desc_class in self.descriptorClasses:
            calculator = desc_class(self)
            temp = calculator(smiles)
            temp.columns = [f"{desc_class.group_name}_{x}" for x in temp.columns]
            self.X = pd.concat([self.X, temp], axis=1)
        return self.X

class QSARModelBuilder(DescriptorBuilderMixIn, CompleteBuilder):

    def getX(self) -> DataFrame:
        return self.X

    def getY(self) -> Series:
        return self.y

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

    def fitAndValidate(
            self,
            X_train : DataFrame,
            y_train : Series,
            X_validated : DataFrame,
            y_validated : Series,
            y_predicted=None,
            perfClass=modelling.models.ModelPerformance,
            *args,
            **kwargs
    ):
        if self.training.mode.name == Algorithm.CLASSIFICATION:
            model = self.algorithmClass(self)
            model.fit(X_train, y_train)
            y_predicted = model.predict(X_validated)
            y_predicted = [1 if x >= 0.5 else 0 for x in y_predicted]
        super().fitAndValidate(X_train, y_train, X_validated, y_validated, y_predicted, perfClass, *args, **kwargs)
