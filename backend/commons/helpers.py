"""
helpers

Created by: Martin Sicho
On: 14-01-20, 17:25
"""

import inspect

def getSubclassesFromModule(base_cls, module):
    ret = []
    for item in inspect.getmembers(module, inspect.isclass):
        if issubclass(item[1], base_cls):
            ret.append(item[1])
    return ret

def findClassInModule(base, module, id_attr, id_attr_val):
    for class_ in base.__subclasses__():
        if hasattr(class_, id_attr):
            if id_attr_val == getattr(class_, id_attr):
                return class_
        else:
            raise Exception("Unspecified ID attribute where required on class: ", repr(class_))