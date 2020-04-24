from django.apps import AppConfig

class GeneratorsConfig(AppConfig):
    name = 'generators'

    def ready(self, force_inspect=False):
        from . import signals
        from modelling import helpers
        helpers.inspectCore("generators", force=force_inspect)

        from commons.helpers import createGroup
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
