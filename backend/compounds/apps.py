from django.apps import AppConfig


class CompoundsConfig(AppConfig):
    name = 'compounds'

    def ready(self):
        from . import signals
        from commons.helpers import createGroup
        from . import models

        createGroup(
            "GenUI_Users",
            [
                models.MolSet,
                models.Activity,
                models.ActivitySet,
                models.ChEMBLActivities,
                models.ChEMBLActivity,
                models.ChEMBLCompounds,
                models.Molecule,
                models.ChEMBLMolecule,
            ]
        )

        createGroup(
            "GenUI_Users",
            [
                models.ActivityTypes,
                models.ActivityUnits,
                models.ChEMBLTarget,
                models.ChEMBLAssay,
                models.MoleculePic,
                models.PictureFormat,

            ],
            permissions=['view']
        )
