import sys

from django.apps import AppConfig
from django.db import transaction


class QsarConfig(AppConfig):
    name = 'qsar'

    def ready(self):
        if len(sys.argv) > 1 and sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate'):
            from .core import bases
            from commons.helpers import getSubclassesFromModule
            from .core import descriptors

            with transaction.atomic():
                for x in getSubclassesFromModule(bases.DescriptorCalculator, descriptors):
                    print(x.getDjangoModel())
