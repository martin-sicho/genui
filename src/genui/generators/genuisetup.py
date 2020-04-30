"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:50 PM
"""

def setup(*args, **kwargs):
    from genui.modelling import helpers
    helpers.inspectCore('genui.generators', force=kwargs['force'])

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
            models.Generator,
            models.GeneratedMolSet,
            models.ModelPerformanceDrugEx,
            models.ModelPerformanceDrugExAgent
        ],
        force=kwargs['force']
    )

