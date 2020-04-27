"""
signals

Created by: Martin Sicho
On: 17-01-20, 10:13
"""
from django.db import IntegrityError
from django.db.models.signals import pre_delete, pre_save
from django.dispatch import receiver

from genui.modelling.models import Model, ModelFile


@receiver(pre_delete, sender=Model, dispatch_uid='on_model_delete_remove_files')
def delete_model_remove_files(sender, instance, using, **kwargs):
    for f in instance.files.all():
        f.file.delete()

@receiver(pre_save, sender=Model, dispatch_uid='on_model_save_checks')
def save_model_checks(sender, instance, using, **kwargs):
    if instance.files.filter(kind=ModelFile.MAIN).count() > 1:
        raise IntegrityError("You can only have one main model file!")

@receiver(pre_save, sender=ModelFile, dispatch_uid="on_file_save")
def save_model_file_callback(sender, instance, using, **kwargs):
    model = Model.objects.get(pk=instance.modelInstance.id)
    model.onFileSave(instance)