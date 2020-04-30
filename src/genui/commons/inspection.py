"""
inspection

Created by: Martin Sicho
On: 4/30/20, 8:55 AM
"""

import importlib
import inspect


def importFromPackage(package, module_name, exception=True):
    package_name = package
    if type(package) != str:
        package_name = package.__name__
    try:
        return importlib.import_module(f'{package_name}.{module_name}')
    except ModuleNotFoundError as exp:
        if exception:
            raise exp
        else:
            return None


def importModuleWithException(module, *args, **kwargs):
    try:
        return importlib.import_module(module, *args, **kwargs)
    except ModuleNotFoundError:
        print(f"WARNING: Failed to find core modelling module: {module}. It will be skipped.")
        return


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


def findSubclassByID(base, module, id_attr : str, id_attr_val : str):
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

    raise Exception(f'Could not find any valid subclass of {base.__name__} in module {module.__name__}.')


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