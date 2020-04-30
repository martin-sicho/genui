"""
inspection

Created by: Martin Sicho
On: 4/30/20, 8:55 AM
"""

import importlib


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
