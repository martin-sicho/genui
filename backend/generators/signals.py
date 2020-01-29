"""
signals

Created by: Martin Sicho
On: 29-01-20, 13:12
"""

from django.db.models.signals import pre_delete
from django.dispatch import receiver

from .models import DrugeExCorpus


@receiver(pre_delete, sender=DrugeExCorpus, dispatch_uid='on_drugexcorpus_delete_remove_files')
def delete_drugexcorpus_remove_files(sender, instance : DrugeExCorpus, using, **kwargs):
    instance.corpusFile.delete()
    instance.vocFile.delete()

