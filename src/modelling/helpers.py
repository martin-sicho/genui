"""
helpers

Created by: Martin Sicho
On: 30-01-20, 13:29
"""
import importlib
import json
import os

from django.db import transaction
from commons.helpers import getSubclassesFromModule, checkInitCondition
from . import models

def importModuleWithExp(module, *args, **kwargs):
    try:
        return importlib.import_module(module, *args, **kwargs)
    except ModuleNotFoundError:
        print(f"WARNING: Failed to find core modelling module: {module}. It will be skipped.")
        return

def inspectCore(referer, core_package="core", modules=("algorithms", "builders", "metrics"), force=False, additional_bases=tuple()):
    if checkInitCondition(force):
        from .core import bases
        base_classes = [bases.Algorithm, bases.ValidationMetric, bases.ModelBuilder] + list(additional_bases)

        with transaction.atomic():
            for module in modules:
                try:
                    module = importlib.import_module(f".{core_package}.{module}", package=referer)
                except ModuleNotFoundError:
                    print(f"Module {referer}.{core_package}.{module} not found. Skipping...")
                    continue

                for base in base_classes:
                    for x in getSubclassesFromModule(base, module):
                        if x == base:
                            continue
                        model = x.getDjangoModel()
                        print(f"Model initialized: {model}")

def createDefaultModels(project, app):
    core_package = importModuleWithExp(f"{app}.core")
    alg_package = importModuleWithExp(f"{app}.core.algorithms")
    models_package = importModuleWithExp(f"{app}.models")

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