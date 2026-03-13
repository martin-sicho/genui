import inspect
import numpy as np
from sklearn import metrics as sklearn_metrics
from qsprpred.models import SklearnMetrics
from qsprpred.models.monitors import BaseMonitor
from genui.utils.inspection import camel_to_snake, snake_to_camel
import traceback
import importlib
from genui.models import models
import re
import gc


def private_state(state):
    return f"_{state}"


class RecordProgressMonitor(BaseMonitor):
    @property
    def validation_index(self):
        return self._validation_index

    @validation_index.setter
    def validation_index(self, value):
        self._validation_index = value

    def __init__(self, progress_callback=None):
        super().__init__()
        self._validation_index = 0
        states = [camel_to_snake(name) for name in dir(self.__class__) if re.match(r"^on[A-Z]", name)]
        for state in states:
            setattr(self, snake_to_camel(state), self._update_callback_method(state))
            setattr(self, private_state(state), False)
            setattr(self.__class__, state, property(
                fget=lambda instance, s=state: getattr(instance, private_state(s)),
                fset=lambda instance, value, s=state: setattr(instance, private_state(s), value)))
        self.progress_callback = progress_callback if progress_callback else lambda: print("Progress updated")

    def _update_callback_method(self, state):
        func = getattr(self, snake_to_camel(state))

        def method(*args, **kwargs):
            if getattr(self, private_state(state)):
                self.progress_callback()
            return func(*args, **kwargs)

        return method

    def onIterationEnd(self, score: float, scores: list[float]):
        super().onIterationEnd(score, scores)
        self.assessments = {}
        gc.collect()

    def onFoldEnd(self, model_fit, fold_predictions):
        super().onFoldEnd(model_fit, fold_predictions)
        self.estimators = {}
        self.fits = {}
        gc.collect()


class MetricsAggregator:
    @property
    def perfClass(self):
        return self._perfClass

    @perfClass.setter
    def perfClass(self, value):
        self._perfClass = value

    def __init__(self, validation_metric, metrics, builder, monitor, perfClass=models.ModelPerformance):
        self.validation_metric = validation_metric
        self.metricFunctions = metrics
        self._perfClass = perfClass
        self.monitor = monitor
        self.builder = builder

    def __call__(self, y_true, y_pred):
        kwargs = {}
        if self.perfClass == models.ModelPerformanceCV:
            kwargs["fold"] = self.monitor.currentFold + 1
            kwargs["validationStrategyIndex"] = self.monitor.validation_index

        for metric in self.metricFunctions[self.monitor.validation_index]:
            if metric != self.validation_metric:
                try:
                    self.builder.saveMetricValue(metric, y_true, y_pred, self.perfClass, **kwargs)
                except Exception as exp:
                    print("Failed to obtain values for metric: ", metric.name)
                    self.builder.errors.append(exp)
                    traceback.print_exc()
        return self.builder.saveMetricValue(self.validation_metric, y_true, y_pred, self.perfClass, **kwargs).value


class CurveMetrics(SklearnMetrics):
    def __init__(self, name):
        super().__init__(self.get_curve_scorer(name))
        self.curve = getattr(sklearn_metrics, name)
        self.score = {'roc_curve': sklearn_metrics.roc_auc_score,
                    'precision_recall_curve': sklearn_metrics.average_precision_score,
                    'det_curve': lambda y_true, y_pred: self._safe_auc_from_curve(self.curve, y_true, y_pred)}

    @staticmethod
    def _safe_auc_from_curve(curve_func, y_true, y_pred):
        x, y = curve_func(y_true, y_pred)[:2]
        if len(x) < 2 or len(y) < 2:
            return np.nan
        return -sklearn_metrics.auc(x, y)

    @staticmethod
    def get_curve_scorer(curve_type):
        def roc_curve(y_true, y_score):
            fpr, tpr, _ = sklearn_metrics.roc_curve(y_true, y_score)
            return sklearn_metrics.auc(fpr, tpr)

        def precision_recall_curve(y_true, y_score):
            precision, recall, _ = sklearn_metrics.precision_recall_curve(y_true, y_score)
            return sklearn_metrics.auc(recall, precision)

        def det_curve(y_true, y_score):
            fpr, fnr, _ = sklearn_metrics.det_curve(y_true, y_score)
            return -sklearn_metrics.auc(fpr, fnr)  # Lower is better, so negate

        scorers = {
            'roc_curve': sklearn_metrics.make_scorer(roc_curve, response_method="predict_proba"),
            'precision_recall_curve': sklearn_metrics.make_scorer(precision_recall_curve,
                                                                  response_method="predict_proba"),
            'det_curve': sklearn_metrics.make_scorer(det_curve, response_method="predict_proba",
                                                     greater_is_better=False)
        }

        if curve_type not in scorers:
            raise ValueError(f"Unsupported curve type: {curve_type}. Choose from {scorers.keys()}.")

        return scorers[curve_type]

    def __call__(self, y_true, y_pred, points=False):
        y_true, y_pred = super().__call__(y_true, y_pred)
        if points:
            return self.curve(y_true, y_pred)
        else:
            return self.score[self.name](y_true, y_pred)

    def _scorerFunc(self, y_true, y_pred):
        return y_true, y_pred


def build_split_instance(obj, dataset=None):
    if isinstance(obj, dict) and "name" in obj:
        # Extract and import the function/class
        full_name = obj["name"]
        module_name, func_name = full_name.rsplit(".", 1)
        module = importlib.import_module(module_name)
        func = getattr(module, func_name)

        # Prepare keyword arguments
        kwargs = {
            k: build_split_instance(v) if isinstance(v, dict) else v
            for k, v in obj.items() if k != "name"
        }
        if dataset is not None and "dataset" in inspect.signature(func).parameters:
            kwargs["dataset"] = dataset

        return func(**kwargs)
    else:
        return obj