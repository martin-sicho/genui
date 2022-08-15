"""
genuisetup

Created by: Martin Sicho
On: 5/3/20, 6:36 PM
"""

PARENT = 'genui.generators'

def setup(*args, **kwargs):
    from genui.utils.init import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.DrugExNet,
            models.DrugExNetValidation,
            models.DrugExNetTraining,
            models.DrugExEnvironment,
            models.DrugExScorer,
            models.DrugExAgent,
            models.DrugExAgentValidation,
            models.DrugExAgentTraining,
            models.DrugEx,
            models.ModelPerformanceDrugEx,
            models.GenUIModelScorer,
            models.PropertyScorer,
            models.ClippedScore
        ],
        force=kwargs['force']
    )

