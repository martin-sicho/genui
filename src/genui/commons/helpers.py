"""
helpers

Created by: Martin Sicho
On: 14-01-20, 17:25
"""
import importlib
import inspect
import logging
import os
import sys

from django.contrib.auth.models import Group, Permission


def getSubclassesFromModule(base_cls, module):
    """
    Fetch all subclasses of a given base class from a module.

    :param base_cls: The base class.
    :param module: The module.
    :return:
    """

    ret = []
    for item in inspect.getmembers(module, inspect.isclass):
        if issubclass(item[1], base_cls):
            ret.append(item[1])
    return ret

def getSubclasses(cls):
    """
    Fetch all existing subclasses of a class.

    :param cls:
    :return:
    """

    return set(cls.__subclasses__()).union([s for c in cls.__subclasses__() for s in getSubclasses(c)])

def findClassbyID(base, module, id_attr : str, id_attr_val : str):
    """
    Function to fetch a given class from a certain module.
    It is identified by both its base class and a value of a
    specific identifying attribute of the class.

    :param base: The base class of the class we are looking for.
    :param module: Module containing the class of interest.
    :param id_attr: The name of the identifying attribute on the class.
    :param id_attr_val: The value of the searched attribute.
    :return:
    """

    all_subclasses = getSubclasses(base)

    for class_ in all_subclasses:
        name = class_.__name__
        if not hasattr(module, name):
            continue

        if not hasattr(class_, id_attr):
            raise Exception("Unspecified ID attribute on a class where it should be defined. Check if the class is properly annotated: ", repr(class_))

        if not (id_attr_val == getattr(class_, id_attr)):
            continue

        return class_

def checkInitCondition(force):
    if 'GENUI_SKIP_INIT' in os.environ and int(os.environ['GENUI_SKIP_INIT']) == 1:
        return False

    return force or (len(sys.argv) > 1 and sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate'))

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

def getFullName(obj, moduleOnly=False):
    module = inspect.getmodule(obj)
    if module is None or module == str.__class__.__module__:
        return obj.__name__ if not moduleOnly else None
    else:
        return module.__name__ + '.' + obj.__name__ if not moduleOnly else module.__name__

def getObjectAndModuleFromFullName(name):
    module = importlib.import_module('.'.join(name.split(".")[0:-1]))
    obj = getattr(module, name.split(".")[-1])
    return obj, module