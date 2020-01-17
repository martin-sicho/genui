"""
signals

Created by: Martin Sicho
On: 17-01-20, 10:13
"""
from django.db.models.signals import pre_delete
from django.dispatch import receiver

from qsar.models import Model


@receiver(pre_delete, sender=Model, dispatch_uid='on_model_delete_remove_files')
def delete_molset_finish_tasks(sender, instance, using, **kwargs):
    instance.modelFile.delete()