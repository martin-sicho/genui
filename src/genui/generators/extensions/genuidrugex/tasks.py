"""
tasks

Created by: Martin Sicho
On: 28-01-20, 13:52
"""
from celery import shared_task

from genui.compounds.models import MolSet, ActivityTypes, Activity
from genui.utils.extensions.tasks.progress import ProgressRecorder
from genui.utils.inspection import getObjectAndModuleFromFullName

from . import models
from .models import DrugExEnvironment, DrugExEnvironmentScores
from .torchutils import cleanup

@shared_task(name="BuildDrugExModel", bind=True, queue='gpu')
def buildDrugExModel(self, model_id, builder_class, model_class):
    # get the builder
    model_class = getattr(models, model_class)
    instance = model_class.objects.get(pk=model_id)
    builder_class = getObjectAndModuleFromFullName(builder_class)[0]
    recorder = ProgressRecorder(self)
    if hasattr(instance, 'parent'):
        builder = builder_class(
            instance,
            instance.parent,
            progress=recorder
        )
    else:
        builder = builder_class(
            instance,
            progress=recorder
        )

    # build the model
    try:
        builder.build()
    except Exception as exp:
        raise exp

    cleanup()
    return {
        "errors" : [repr(x) for x in builder.errors],
        "DrExModelName" : instance.name,
        "DrExModelID" : instance.id,
    }

@shared_task(name="calculateEnvironment", bind=True)
def calculateEnvironment(self, environment_id, molsets_ids, use_modifiers):
    environment = DrugExEnvironment.objects.get(pk=environment_id)
    instance = environment.getInstance(use_modifiers=use_modifiers)
    activity_sets = []
    for molset_id in molsets_ids:
        molset = MolSet.objects.get(pk=molset_id)
        molecules = molset.molecules.all()
        scores = instance(molset.allSmiles)
        assert len(scores) == len(molecules)
        activity_set = DrugExEnvironmentScores(
            molecules=molset,
            environment=environment,
            project=environment.project,
            name=environment.name,
            description=f"Activity set created from environment: {environment.name} for compound set: {molset.name}."
        )
        activity_set.save()
        activity_sets.append(activity_set.id)

        drops = ['VALID']
        if not use_modifiers:
            drops.append('DESIRE')
        scores.drop(drops, axis=1, inplace=True)
        columns = scores.columns
        for idx, row in scores.iterrows():
            molecule = molecules[idx]
            for col in columns:
                value = row[col]
                if value:
                    atype = ActivityTypes.objects.get_or_create(value=f"{environment.name}_{col}")[0]

                    activity = Activity(
                        value=value,
                        type=atype,
                        source=activity_set,
                        molecule=molecule
                    )
                    activity.save()

    return {
        "environment" : environment.id,
        "activitySets" : activity_sets,
    }
