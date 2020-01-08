"""
signals

Created by: Martin Sicho
On: 1/5/20, 5:40 PM
"""
from django.db.models.signals import pre_delete
from django.dispatch import receiver

from .models import MolSet

@receiver(pre_delete, sender=MolSet)
def delete_molset(sender, instance, using):
    # TODO: find all running tasks connected tot his instance and revoke them
    print('About to delete: ', repr(instance))
