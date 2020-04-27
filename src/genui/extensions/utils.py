"""
utils

Created by: Martin Sicho
On: 4/27/20, 6:35 PM
"""
import importlib


def discover_extensions(packages):
    ret = []
    for package_name in packages:
        try:
            package = importlib.import_module(package_name)
        except ModuleNotFoundError as exp:
            print(f"WARNING: Failed to find extension module: {package_name}.")
            raise exp

        subs = package.__all__
        for sub in subs:
            ret.append(f"{package_name}.{sub}")
    return ret

