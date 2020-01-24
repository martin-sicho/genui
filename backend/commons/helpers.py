"""
helpers

Created by: Martin Sicho
On: 14-01-20, 17:25
"""

import inspect

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

def findClassInModule(base, module, id_attr : str, id_attr_val : str):
    """
    Function to fetch a given class from a certain module.
    It is identified by both its base class and a value of a
    specific identifying attribute of the class.

    :param base: The base class of the class we are looking for.
    :param module: This argument is not used in the function, but it ensures that the module where the class is defined will be imported before we attempt to search for the class.
    :param id_attr: The name of the identifying attribute on the class.
    :param id_attr_val: The value of the searched attribute.
    :return:
    """

    for class_ in base.__subclasses__():
        if hasattr(class_, id_attr):
            if id_attr_val == getattr(class_, id_attr):
                return class_
        else:
            raise Exception("Unspecified ID attribute on a class where it should be defined. Check if the class is properly annotated: ", repr(class_))