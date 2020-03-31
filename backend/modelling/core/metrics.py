"""
metrics

Created by: Martin Sicho
On: 24-01-20, 15:04
"""
from pandas import Series
from sklearn import metrics

from . import bases


class MCC(bases.ValidationMetric):
    name = "MCC"
    description = "Matthew's Correlation Coefficient"
    modes = [bases.Algorithm.CLASSIFICATION]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        return metrics.matthews_corrcoef(true_vals, predicted_vals)