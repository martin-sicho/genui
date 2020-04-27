from django.apps import AppConfig


class CompoundsConfig(AppConfig):
    name = 'genui.compounds'

    def ready(self, force=False):
        from . import signals
        from genui.commons.helpers import createGroup
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
            ],
            force=force
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
            permissions=['view'],
            force=force
        )
