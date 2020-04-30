"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:43 PM
"""

def setup(*args, **kwargs):
    from . import signals
    from . import helpers
    helpers.inspectCore('genui.modelling', force=kwargs['force'])

    from genui.utils.init import createGroup
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
        force=kwargs['force']
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
        permissions=['view'],
        force=kwargs['force']
    )

