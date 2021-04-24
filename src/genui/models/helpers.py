"""
helpers

Created by: Martin Sicho
On: 30-01-20, 13:29
"""
import importlib
import json
import os

from django.db import transaction
from genui.utils.init import checkInitCondition
from . import models
from genui.utils.inspection import importModuleWithException, getSubclassesFromModule


def discoverGenuiModels(container, core_package="genuimodels", modules=("algorithms", "builders", "metrics"), force=False, additional_bases=tuple()):
    if checkInitCondition(force):
        from .genuimodels import bases
        base_classes = [bases.Algorithm, bases.ValidationMetric, bases.ModelBuilder] + list(additional_bases)

        with transaction.atomic():
            for module in modules:
                try:
                    module = importlib.import_module(f"{container}.{core_package}.{module}")
                except ModuleNotFoundError as err:
                    # print(f"Module {container}.{core_package}.{module} failed to import. It will be skipped. Reason: {err}")
                    if f"{container}.{core_package}" not in repr(err):
                        raise err

                for base in base_classes:
                    for x in getSubclassesFromModule(base, module):
                        if x == base:
                            continue
                        model = x.getDjangoModel(corePackage=f"{container}.{core_package}")
                        print(f"Django model instance initialized for '{model}' from module: '{module.__name__}'")

def createDefaultModels(project, app):
    core_package = importModuleWithException(f"{app}.genuimodels", message=False, throw=True)
    alg_package = importModuleWithException(f"{app}.genuimodels.algorithms", message=False, throw=True)
    models_package = importModuleWithException(f"{app}.models", message=False, throw=True)

    if not (core_package and alg_package and models_package):
        return

    path = os.path.dirname(core_package.__file__)
    path = os.path.join(path, "default_models")

    if not os.path.exists(path):
        return

    ret = []
    for subdir in os.listdir(path):
        model_dir = os.path.join(path, subdir)
        metadata = json.load(open(os.path.join(model_dir, "model.json")))
        model_class = getattr(models_package, metadata["modelClass"])
        if not model_class.objects.filter(name=metadata["name"], project=project).exists():
            instance = model_class.objects.create(
                name=metadata["name"],
                description=metadata["description"],
                project=project,
                builder=models.ModelBuilder.objects.get(name=metadata["builderClass"])
            )

            ts_data = metadata["trainingStrategy"]
            training_strat_class = getattr(models_package, ts_data["className"])
            alg_class = getattr(alg_package, ts_data["algorithmClass"])
            ts = training_strat_class.objects.create(
                modelInstance=instance,
                algorithm=models.Algorithm.objects.get(name=alg_class.__name__),
                mode=models.AlgorithmMode.objects.get(name=ts_data['mode']),
            )

            if ts.mode.name == "generator":
                gen_class = getattr(models_package, metadata["generatorClass"])
                gen_class.objects.create(
                    agent=instance,
                    name=instance.name,
                    description=instance.description,
                    project=instance.project
                )

            for file_ in metadata["files"]:
                file_path = os.path.join(model_dir, file_["path"])
                models.ModelFile.create(
                    instance
                    , os.path.basename(file_path)
                    , open(file_path, mode="rb")
                    , kind=file_["kind"]
                    , note=file_["note"] if "note" in file_ else None
                )

            ret.append(instance)

    return ret