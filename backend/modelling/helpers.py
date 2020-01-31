"""
helpers

Created by: Martin Sicho
On: 30-01-20, 13:29
"""
import importlib

from django.db import transaction
import sys
from commons.helpers import getSubclassesFromModule

def inspectCore(referer, core_package="core", modules=("algorithms", "builders", "metrics"), force=False, additional_bases=tuple()):
    if force or (len(sys.argv) > 1 and sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate', "test")):
        from .core import bases
        base_classes = [bases.Algorithm, bases.ValidationMetric] + list(additional_bases)

        with transaction.atomic():
            for module in modules:
                try:
                    module = importlib.import_module(f".{core_package}.{module}", package=referer)
                except ModuleNotFoundError:
                    print(f"Module {referer}.{core_package}.{module} not found. Skipping...")
                    continue

                for base in base_classes:
                    for x in getSubclassesFromModule(base, module):
                        if x == base:
                            continue
                        model = x.getDjangoModel()
                        print(f"Model initialized: {model}")

