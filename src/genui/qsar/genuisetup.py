
def setup(*args, **kwargs):
    from genui.models import helpers
    from .genuimodels import bases
    from genui import apps
    for app in apps.all_():
        helpers.discoverGenuiModels(app, force=kwargs['force'], modules=["embeddings"], additional_bases=[bases.EmbeddingCalculator])

    from genui.utils.init import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.QSARModel,
            models.ModelActivity,
            models.ModelActivitySet,
            models.QSARTrainingStrategy,
        ],
        force=kwargs['force']
    )

    createGroup(
        "GenUI_Users",
        [
            models.EmbeddingCalculator,
        ],
        permissions=['view'],
        force=kwargs['force']
    )
