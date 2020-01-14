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