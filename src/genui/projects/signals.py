"""
signals

Created by: Martin Sicho
On: 3/30/20, 12:16 PM
"""
from django.contrib.auth.models import User, Group
from django.db.models.signals import post_save
from django.dispatch import receiver


@receiver(post_save, sender=User, dispatch_uid='on_user_save_add_user_to_users_group')
def add_user_to_users_group(sender, instance, created, **kwargs):
    if created:
        users_group = Group.objects.get_or_create(name="GenUI_Users")[0]
        instance.groups.add(users_group)
        # print(f"New User '{instance.username}' added to '{users_group.name}'.")

