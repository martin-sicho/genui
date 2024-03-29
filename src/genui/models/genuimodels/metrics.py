"""
metrics

Created by: Martin Sicho
On: 24-01-20, 15:04
"""
from pandas import Series
from sklearn import metrics
from sklearn.metrics import r2_score, mean_squared_error, roc_curve, auc

from . import bases
from .. import models

class ROC(bases.ValidationMetric):
    name = "ROC"
    description = "The area under ROC curve."
    modes = [bases.Algorithm.CLASSIFICATION]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        fpr, tpr, thrs = self.getCurve(true_vals, predicted_vals)
        return auc(fpr, tpr)

    def getCurve(self, true_vals, scores):
        return roc_curve(true_vals, scores)

    def save(
            self,
            true_vals : Series,
            predicted_vals : Series,
            perfClass=models.ModelPerformance,
            **kwargs
    ):
        roc_auc = super().save(true_vals, predicted_vals, perfClass, **kwargs)
        for fpr, tpr, trh in zip(*self.getCurve(true_vals, predicted_vals)):
            models.ROCCurvePoint.objects.create(
                metric=models.ModelPerformanceMetric.objects.get(name=self.name),
                model=self.builder.instance,
                value=tpr,
                fpr=fpr,
                auc=roc_auc
            )
        return roc_auc


class MCC(bases.ValidationMetric):
    name = "MCC"
    description = "Matthew's Correlation Coefficient"
    modes = [bases.Algorithm.CLASSIFICATION]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        return metrics.matthews_corrcoef(true_vals, self.probasToClasses(predicted_vals))

class R2(bases.ValidationMetric):
    name = "R2"
    description = "R^2 (coefficient of determination) regression score function. As implemented in scikit-learn: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html#sklearn.metrics.r2_score."
    modes = [bases.Algorithm.REGRESSION]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        return r2_score(true_vals, predicted_vals)

class MSE(bases.ValidationMetric):
    name = "MSE"
    description = "Mean squared error regression loss. As implemented in scikit-learn: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.mean_squared_error.html#sklearn.metrics.mean_squared_error."
    modes = [bases.Algorithm.REGRESSION]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        return mean_squared_error(true_vals, predicted_vals)