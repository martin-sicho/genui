from django.apps import AppConfig

class ModellingConfig(AppConfig):
    name = 'genui.modelling'

    def ready(self, force_inspect=False):
        from . import signals
        from . import helpers
        helpers.inspectCore('genui.modelling', force=force_inspect)

        from genui.commons.helpers import createGroup
        from . import models

        createGroup(
            "GenUI_Users",
            [
                models.Model,
                models.TrainingStrategy,
                models.ValidationStrategy,
                models.BasicValidationStrategy,
                models.ModelFile,
                models.ModelParameterValue,
                models.ModelParameterBool,
                models.ModelParameterFloat,
                models.ModelParameterInt,
                models.ModelParameterStr,
                models.ModelPerformance,
                models.ModelPerfomanceNN,
                models.ModelPerformanceCV,
            ],
            force=force_inspect
        )

        createGroup(
            "GenUI_Users",
            [
                models.Algorithm,
                models.AlgorithmMode,
                models.ModelBuilder,
                models.ModelFileFormat,
                models.ModelParameter,
                models.ModelPerformanceMetric,
            ],
            permissions=['view']
        )
