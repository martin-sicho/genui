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
            models.DrugEx,
            models.DrugExNet,
            models.DrugExNetTrainingStrategy,
            models.DrugExValidationStrategy,
            models.DrugExAgent,
            models.DrugExAgentTrainingStrategy,
            models.DrugExAgentValidationStrategy,
            models.ModelPerformanceDrugEx,
            models.ModelPerformanceDrugExAgent
        ],
        force=kwargs['force']
    )

