from django.apps import AppConfig

class GeneratorsConfig(AppConfig):
    name = 'genui.generators'

    def ready(self, force_inspect=False):
        from . import signals
        from genui.modelling import helpers
        helpers.inspectCore('genui.generators', force=force_inspect)

        from genui.commons.helpers import createGroup
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
            force=force_inspect
        )
