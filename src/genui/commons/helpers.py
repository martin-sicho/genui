"""
helpers

Created by: Martin Sicho
On: 14-01-20, 17:25
"""
import logging

from django.contrib.auth.models import Group, Permission

from genui.utils.init import checkInitCondition


def createGroup(
        groupName
        , models
        , permissions=('add', 'change', 'delete', 'view')
        , overwrite=True
        , appendPermissions=True
        , force=False
):
    if checkInitCondition(force):
        group, created = Group.objects.get_or_create(name=groupName)
        if not overwrite and not created:
            raise Exception(f"Group {groupName} already exists and overwrite is off.")
        if not appendPermissions:
            group.permissions.clear()

        for model in models:
            for permission in permissions:
                codename = f"{permission}_{model.__name__.lower()}"
                print(f"Creating permission for group {groupName}: {codename}")

                try:
                    model_add_perm = Permission.objects.get(codename=codename)
                except Permission.DoesNotExist:
                    logging.warning("Permission not found with name '{}'.".format(codename))
                    continue

                group.permissions.add(model_add_perm)

