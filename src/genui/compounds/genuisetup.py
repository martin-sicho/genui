"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:38 PM
"""
import importlib

from genui.utils.inspection import getSubclassesFromModule


def discoverExporters(app, exporters_module="exporters"):
    from .exporters.base import BaseMolSetExporter
    from .models import MolSetExporter

    name = f"{app}.{exporters_module}"
    module = None
    try:
        module = importlib.import_module(name)
    except ModuleNotFoundError as err:
        if name not in repr(err):
            raise err
        else:
            return

    for exporter in getSubclassesFromModule(BaseMolSetExporter, module):
        if exporter == BaseMolSetExporter:
            continue
        else:
            instance = MolSetExporter.objects.get_or_create(
                name=exporter.name
            )[0]
            instance.classPath = f"{name}.{exporter.__name__}"
            instance.save()
            print(f"Found molecule set exporter: {instance}")


def setup(*args, **kwargs):
    from genui.utils.init import createGroup
    from genui import apps
    from . import models
    from . import signals

    for app in apps.all_():
        discoverExporters(app)

    createGroup(
        "GenUI_Users",
        [
            models.MolSet,
            models.Activity,
            models.ActivitySet,
            models.Molecule,
            models.MolSetFile,
            models.MolSetExport
        ],
        force=kwargs['force']
    )

    createGroup(
        "GenUI_Users",
        [
            models.ActivityTypes,
            models.ActivityUnits,
            models.MoleculePic,
            models.PictureFormat,
            models.MolSetExporter
        ],
        permissions=['view'],
        force=kwargs['force']
    )

