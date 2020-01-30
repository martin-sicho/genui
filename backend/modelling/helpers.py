"""
helpers

Created by: Martin Sicho
On: 30-01-20, 13:29
"""
import importlib

from django.db import transaction
import sys
from commons.helpers import getSubclassesFromModule

def inspectCore(referer, core_package="core", modules=("algorithms", "builders", "metrics"), force=False):
    if force or (len(sys.argv) > 1 and sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate', "test")):
        from .core import bases

        with transaction.atomic():
            for module in modules:
                try:
                    module = importlib.import_module(f".{core_package}.{module}", package=referer)
                except ModuleNotFoundError:
                    print(f"Module {referer}.{core_package}.{module} not found. Skipping...")
                    continue
                for x in getSubclassesFromModule(bases.Algorithm, module):
                    if x == bases.Algorithm:
                        continue
                    print(f"Model initialized: {x.getDjangoModel()}")
                    print(f"Found parameters for {x.name}: {x.getParams()}")

                for x in getSubclassesFromModule(bases.ValidationMetric, module):
                    if x == bases.ValidationMetric:
                        continue
                    print(f"Model initialized: {x.getDjangoModel()}")

