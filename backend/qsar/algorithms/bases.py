"""
algorithms

Created by: Martin Sicho
On: 14-01-20, 10:16
"""
import uuid
from django.core.files.base import ContentFile

import modelling.models
from commons.helpers import findClassInModule
from modelling.algorithms.bases import Algorithm, ValidationMetric
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


class QSARModelBuilder:
    @staticmethod
    def findAlgorithmClass(name):
        from . import algorithms
        return findClassInModule(Algorithm, algorithms, "name", name)

    @staticmethod
    def findMetricClass(name):
        from . import metrics
        return findClassInModule(ValidationMetric, metrics, "name", name)

    @staticmethod
    def findDescriptorClass(name):
        from . import descriptors
        return findClassInModule(DescriptorCalculator, descriptors, "group_name", name)

    def __init__(
            self,
            instance : models.QSARModel,
            progress = None,
            onFitCall=None
    ):
        self.instance = instance
        self.training = self.instance.trainingStrategy
        self.validation = self.instance.validationStrategy
        self.descriptorClasses = [self.findDescriptorClass(x.name) for x in self.training.descriptors.all()]
        self.algorithmClass = self.findAlgorithmClass(self.training.algorithm.name)
        self.metricClasses = [self.findMetricClass(x.name) for x in self.validation.metrics.all()]
        self.model = None
        self.onFitCall = onFitCall
        self.molset = self.instance.molset
        self.progress = progress
        self.X = None
        self.y = None
        self.progressStages = []
        self.currentProgress = 0
        self.errors = []

    def recordProgress(self):
        if self.currentProgress < len(self.progressStages):
            if self.progress:
                self.progress.set_progress(
                    self.currentProgress
                    , len(self.progressStages)
                    , description=self.progressStages[self.currentProgress]
                )
            print(self.progressStages[self.currentProgress])
        else:
            self.errors.append(Exception("Incorrect progress count detected."))
        self.currentProgress += 1
        print(f"{self.currentProgress}/{len(self.progressStages)}")

    def calculateDescriptors(self, mols):
        smiles = [x.canonicalSMILES for x in mols]
        if not self.X:
            self.X = DataFrame()
            for desc_class in self.descriptorClasses:
                calculator = desc_class(self)
                temp = calculator(smiles)
                temp.columns = [f"{desc_class.group_name}_{x}" for x in temp.columns]
                self.X = pd.concat([self.X, temp], axis=1)
            return self.X

    def saveActivities(self):
        if not self.y:
            activity_set = self.instance.molset.activities.get()
            compounds, activities = activity_set.cleanForModelling()
            activities = Series(activities)
            if self.training.mode.name == Algorithm.CLASSIFICATION:
                activity_thrs = self.training.activityThreshold
                activities = activities.apply(lambda x : 1 if x >= activity_thrs else 0)
            self.y = activities
            return self.y, compounds

    def fit(self, callback=None) -> models.QSARModel:
        # TODO: sanity check if length of X and y are the same
        self.model = self.algorithmClass(self, callback if callback else self.onFitCall)
        self.model.fit(self.X, self.y)
        self.saveFile()
        return self.instance

    def fitValidate(self) -> models.QSARModel:
        ret = self.fit()
        self.validate(self.model, self.X, self.y)
        return ret

    def validate(self, model, X : DataFrame, y_truth : Series, perfClass=modelling.models.ModelPerformance, **kwargs):
        predictions = model.predict(X)[:,1]
        if self.training.mode.name == Algorithm.CLASSIFICATION:
            predictions = [1 if x >= 0.5 else 0 for x in predictions]
        for metric_class in self.metricClasses:
            performance = metric_class(self)(y_truth, predictions)
            perfClass.objects.create(
                    metric=modelling.models.ModelPerformanceMetric.objects.get(name=metric_class.name),
                    value=performance,
                    model=self.instance,
                    **kwargs
                )

    def saveFile(self):
        name = f"{self.algorithmClass.name}{self.instance.id}_project{self.instance.project.id}_{uuid.uuid1()}"
        extension = self.model.fileFormat.fileExtension
        self.instance.modelFile.save(name + extension, ContentFile('Dummy file for {0}'.format(name)))
        path = self.instance.modelFile.path.replace(extension, '')
        self.model.serialize(path)