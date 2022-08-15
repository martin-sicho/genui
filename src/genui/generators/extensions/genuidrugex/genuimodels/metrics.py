"""
metrics

Created by: Martin Sicho
On: 27-01-20, 11:08
"""
from pandas import Series

from genui.models.genuimodels.bases import ValidationMetric, Algorithm


class SMILESErrorRate(ValidationMetric):
    name = "SMILES_ER"
    description = "Percentage of invalid smiles in the generated structures."
    modes = [Algorithm.GENERATOR]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        # TODO: implement this
        raise NotImplementedError("This is not yet supported.")

class SMILESUniqueRate(ValidationMetric):
    name = "SMILES_UQR"
    description = "Percentage of valid unique smiles that scored above the decision threshold in the predicted activity values."
    modes = [Algorithm.GENERATOR]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        # TODO: implement this
        raise NotImplementedError("This is not yet supported.")

class MeanDrExDesirability(ValidationMetric):
    name = "DrExDesire"
    description = "Ratio of the desired molecules generated in the set."
    modes = [Algorithm.GENERATOR]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        raise NotImplementedError("This method is not intended to be called, because DrugEx has this calculation built in.")

class DrugExLoss(ValidationMetric):
    name = "DrExLoss"
    description = "Value of the DrugEx loss function."
    modes = [Algorithm.GENERATOR]

    def __call__(self, true_vals: Series, predicted_vals: Series):
        raise NotImplementedError("This method is not intended to be called, because DrugEx has this calculation built in.")