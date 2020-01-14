"""
metrics

Created by: Martin Sicho
On: 14-01-20, 14:06
"""

from pandas import Series
from sklearn import metrics

from qsar import models
from . import bases


class MCC(bases.ValidationMetric):
    name = "MCC"
    description = ""

    @staticmethod
    def getDjangoModel():
        ret = models.ModelPerformanceMetric.objects.get_or_create(
            name=MCC.name
        )[0]
        ret.description = MCC.description
        ret.save()

    def __call__(self, true_vals: Series, predicted_vals: Series):
        return metrics.matthews_corrcoef(true_vals, predicted_vals)