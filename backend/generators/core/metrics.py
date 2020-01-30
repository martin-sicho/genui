"""
metrics

Created by: Martin Sicho
On: 27-01-20, 11:08
"""
from pandas import Series

from modelling.core.bases import ValidationMetric


class SMILESErrorRate(ValidationMetric):
    name = "SMILES_ER"
    description = "Percentage of invalid smiles in the predicted values"

    def __call__(self, true_vals: Series, predicted_vals: Series):
        # TODO: implement this
        raise NotImplementedError("This is not yet supported.")

class DrugExLoss(ValidationMetric):
    name = "DrExLoss"
    description = "Value of the DrugEx loss function."

    def __call__(self, true_vals: Series, predicted_vals: Series):
        raise NotImplementedError("This method is not intended to be called, because DrugEx has this calculation built in.")