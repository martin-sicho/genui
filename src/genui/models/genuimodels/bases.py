import logging
import traceback
import weakref
from abc import ABC, abstractmethod

import joblib
import numpy as np

from django.db import transaction
from pandas import DataFrame, Series
from sklearn import metrics

from django.core.files.base import ContentFile
from qsprpred.data import QSPRDataset
from qsprpred.models.assessment.metrics import classification, regression, scikit_learn

from genui.models import models
from genui.utils.exceptions import GenUIException
from genui.utils.inspection import findSubclassByID, importFromPackage
from genui.models.models import ModelFile
from genui.qsar.genuimodels.qsprpred_utils import CurveMetrics


class Algorithm(ABC):
    # TODO: use a metaclass to initialize these (the classmethods below have to be called explicitly at the moment)
    name = None
    parameters = {}
    CLASSIFICATION = 'classification'
    MULTICLASS = 'multiclass'
    REGRESSION = 'regression'
    GENERATOR = 'generator'
    MAP = 'map'

    django_model = None
    django_parameters = []
    django_file_formats = []
    django_modes = []

    @classmethod
    def attachModesToModel(cls, model, modes):
        model.validModes.clear()
        for mode in modes:
            mode = models.AlgorithmMode.objects.get_or_create(name=mode)[0]
            model.validModes.add(mode)
        model.save()
        return model.validModes

    @classmethod
    def getModes(cls):
        return [cls.CLASSIFICATION, cls.REGRESSION]

    @staticmethod
    def attachToInstance(instance, items, field):
        field.clear()
        for item in items:
            field.add(item)
        instance.save()

    @classmethod
    def getFileFormats(cls, attach_to=None):
        formats = [models.ModelFileFormat.objects.get_or_create(
            fileExtension=".joblib.gz",
            description="A compressed joblib file."
        )[0], ]
        if attach_to:
            cls.attachToInstance(attach_to, formats, attach_to.fileFormats)
        return formats

    @classmethod
    def getDjangoModel(cls, corePackage=None, update=False):
        # TODO: this should go to the init of the metaclass
        if not cls.name:
            print('This class has invalid name attribute. No django model can be provided for: ', cls.__name__)
            return

        ret, ret_created = models.Algorithm.objects.get_or_create(name=cls.name)

        # just return if we are not setting up a new instance
        if not ret_created and not update:
            return ret

        if corePackage:
            ret.corePackage = corePackage
            ret.save()

        cls.django_modes = cls.attachModesToModel(ret,
                                                  cls.getModes())  # TODO: this should use the same pattern as the file formats method
        cls.django_file_formats = cls.getFileFormats(attach_to=ret)
        cls.django_model = ret
        cls.django_parameters = cls.getParams()
        return ret

    @classmethod
    def getParams(cls):
        ret = []
        current_params = models.ModelParameter.objects.filter(
            algorithm=cls.django_model
        )
        missing_params = {x for x in models.ModelParameter.objects.filter(
            algorithm=cls.django_model
        ).all()}
        for param_name in cls.parameters:
            param_type = cls.parameters[param_name]["type"]
            with transaction.atomic():
                param, created = current_params.get_or_create(
                    name=param_name,
                    contentType=param_type,
                    algorithm=cls.django_model
                )

                default_value = cls.parameters[param_name]["defaultValue"]
                if created:
                    print(f"Creating default value: {cls.__name__}.{param_name} = {default_value}")
                    param.defaultValue = models.PARAM_VALUE_CTYPE_TO_MODEL_MAP[param_type].objects.create(
                        parameter=param,
                        value=default_value
                    )
                    param.save()
                else:
                    if param.defaultValue.value != default_value:
                        print(f'Changing default value: {param_name}={default_value}')
                        param.defaultValue.value = default_value
                        param.defaultValue.save()

                ret.append(param)
                if param in missing_params:
                    missing_params.remove(param)

        if missing_params:
            for param in missing_params:
                logging.warning(
                    f"Parameter {param} no longer present for algorithm {cls}. It will be removed from the database and from all prior models that use it.")
                param.delete()

        return ret

    def __init__(self, builder, callback=None):
        self._builder = weakref.ref(builder)
        self.instance = self.builder.instance
        self.trainingInfo = self.builder.training
        self.validationInfo = self.builder.validations
        self.params = {x.parameter.name: x.value for x in self.trainingInfo.parameters.all()}
        self.mode = self.trainingInfo.mode
        self.callback = callback
        self._model = None

    @property
    def builder(self):
        ret = self._builder()
        if ret:
            return ret
        else:
            raise LookupError("Builder was destroyed before being referenced!")

    def getSerializer(self):
        return lambda filename: joblib.dump(
            self.model
            , filename
        )

    def serialize(self, filename):
        self.getSerializer()(filename)

    def getDeserializer(self):
        return lambda filename: joblib.load(filename)

    def deserialize(self, filename):
        self._model = self.getDeserializer()(filename)
        return self

    @property
    @abstractmethod
    def model(self):
        pass

    @abstractmethod
    def fit(self, X: DataFrame, y: Series):
        pass

    @abstractmethod
    def predict(self, X: DataFrame) -> Series:
        pass


class ModelBuilder(ABC):

    @classmethod
    def getDjangoModel(cls, corePackage=None, update=False):
        ret, ret_created = models.ModelBuilder.objects.get_or_create(
            name=cls.__name__
        )

        # just return if we are not setting up a new instance
        if not ret_created and not update:
            return ret

        if corePackage:
            ret.corePackage = corePackage
            ret.save()
        return ret

    def findAlgorithmClass(self, name, corePackage=None):
        if not corePackage:
            corePackage = self.corePackage
        return findSubclassByID(
            Algorithm,
            importFromPackage(corePackage, "algorithms"),
            "name",
            name
        )

    def getMetricFunction(self, name):
        if "curve" in name:
            return CurveMetrics(name)
        else:
            return scikit_learn.SklearnMetrics(name)

    def __init__(
            self,
            instance: models.Model,
            progress=None,
            onFit=None
    ):
        self.instance = instance

        self.training = self.instance.trainingStrategy
        self.algorithmPackage = self.training.algorithm.corePackage
        self.algorithmClass = self.findAlgorithmClass(
            self.training.algorithm.name,
            self.algorithmPackage if self.algorithmPackage else None
        )
        self.onFit = onFit

        self.validations = self.training.validationStrategies.all()
        self.hyper_param_opt = self.training.hyperParamOptStrategies.first()
        self.metricFunctions = []
        for validation in self.validations:
            current_metrics = []
            for metric in validation.metrics:
                current_metrics.append(self.getMetricFunction(metric))
            self.metricFunctions.append(current_metrics)
        if self.hyper_param_opt:
            self.hyper_param_aggregator = getattr(np, self.hyper_param_opt.scoreAggregation)
            self.hypo_metric = self.getMetricFunction(self.hyper_param_opt.metric)
            self.search_space = self.prepare_seachSpace(self.hyper_param_opt)
        self.progress = progress
        self.errors = []

        self._model = None

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

    @property
    def model(self) -> Algorithm:
        if self._model is None:
            self._model = self.algorithmClass(self, self.onFit)
            if self.instance.modelFile:
                self._model.deserialize(self.instance.modelFile.path)
        return self._model

    @model.setter
    def model(self, val):
        self._model = val

    @abstractmethod
    def build(self) -> models.Model:
        pass

    @abstractmethod
    def validate(self, validation_strategy):
        # Implement validation logic here
        pass

    def saveFile(self):
        if not self.instance.modelFile:
            model_format = self.training.algorithm.fileFormats.all()[
                0]  # FIXME: this should be changed once we expose the file formats in the training strategy
            ModelFile.create(
                self.instance,
                f'main.{model_format.fileExtension}',
                ContentFile('placeholder'),
                kind=ModelFile.MAIN,
                note=f'{self.training.algorithm.name}_main'
            )
        path = self.instance.modelFile.path
        self.model.serialize(path)

    def prepare_seachSpace(self, hyper_param_opt):
        orig_search_space = hyper_param_opt.searchSpace
        search_space = {}
        for param in orig_search_space:
            param_name = param['name']
            param_type = param['type']
            param_value = param['value']
            if param_type == 'range':
                param_value = range(*param_value)
            elif param_type == 'int' or param_type == 'float':
                param_value = [param_type, *param_value]
            elif param_type == 'categorical':
                param_value = [param_type, param_value]
            search_space[param_name] = param_value
        return search_space


class ProgressMixIn:

    def __init__(self, instance, progress, *args, **kwargs):
        super().__init__(instance, progress, *args, **kwargs)

        self.progress = progress
        self.progressStages = []
        self.currentProgress = 0
        self.errors = self.errors if hasattr(self, "errors") else []

    def recordProgress(self):
        if self.currentProgress < len(self.progressStages):
            if self.progress:
                self.progress.set_progress(
                    self.currentProgress + 1
                    , len(self.progressStages)
                    , description=self.progressStages[self.currentProgress]
                )
            print(self.progressStages[self.currentProgress])
        else:
            self.errors.append(Exception("Incorrect progress count detected."))
        self.currentProgress += 1
        print(f"{self.currentProgress}/{len(self.progressStages)}")


class ValidationMixIn:

    def fitAndValidate(
            self,
            X_train: DataFrame,
            y_train: Series,
            X_validated: DataFrame,
            y_validated: Series,
            y_predicted=None,
            perfClass=models.ModelPerformance,
            *args,
            **kwargs
    ):
        if not y_predicted:
            model = self.algorithmClass(self)
            model.fit(X_train, y_train)
            y_predicted = model.predict(X_validated)
        self.validate(y_validated, y_predicted, perfClass, *args, **kwargs)

    def _saveMetricValue(self, metric, y_true, y_predicted, perfClass=models.ModelPerformance, *args, **kwargs):
        value = metric(y_true, y_predicted)
        return perfClass.objects.create(
            metric=metric.name,
            value=value,
            model=self.instance,
            **kwargs
        )

    def saveMetricValue(self, metric, y_true, y_predicted, perfClass=models.ModelPerformance, *args, **kwargs):
        if "curve" in metric.name:
            return self.saveCurvePoints(metric, y_true, y_predicted, perfClass, *args, **kwargs)
        else:
            return self._saveMetricValue(metric, y_true, y_predicted, perfClass, *args, **kwargs)

    def saveCurvePoints(self, metric, y_true, y_predicted, perfClass=models.ModelPerformance, *args, **kwargs):
        dependent, independent, _ = metric(y_true, y_predicted, True)
        perf_object = self._saveMetricValue(metric, y_true, y_predicted, perfClass, *args, **kwargs)
        for ind, dep in zip(independent, dependent):
            models.MetricCurvePoint.objects.create(
                metric=metric.name,
                model=self.instance,
                value=dep,
                independent=ind,
                auc=perf_object,
            )
        return perf_object

    def validate(self, y_validated, y_predicted, perfClass=models.ModelPerformance, *args, **kwargs):
        metric_functions = set(self.metricFunctions) if not isinstance(self.metricFunctions[0], list) \
            else set([mc for mcs in self.metricFunctions for mc in mcs])
        for metric in metric_functions:
            try:
                self.saveMetricValue(metric, y_validated, y_predicted, *args, **kwargs)
            except Exception as exp:
                print("Failed to obtain values for metric: ", metric.name)
                self.errors.append(exp)
                traceback.print_exc()


class PredictionMixIn:

    def predict(self, dataset: DataFrame | QSPRDataset | np.ndarray = None) -> Series:
        if self.model:
            if dataset is None:
                return self.predictMols(dataset)
            return self.model.predict(dataset)
        else:
            raise ModelNotFittedException("The model is not trained or loaded. Invalid call to 'predict'.")


class ModelNotFittedException(GenUIException):
    pass
