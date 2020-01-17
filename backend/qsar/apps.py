import sys

from django.apps import AppConfig
from django.db import transaction


class QsarConfig(AppConfig):
    name = 'qsar'

    def ready(self):
        from . import signals
        if sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate'):
            from .algorithms import bases
            from commons.helpers import getSubclassesFromModule
            from .algorithms import algorithms as algs
            from .algorithms import metrics
            from .algorithms import descriptors

            with transaction.atomic():
                for x in getSubclassesFromModule(bases.Algorithm, algs):
                    print(x.getParams())
                    print(x.getDjangoModel())

                for x in getSubclassesFromModule(bases.ValidationMetric, metrics):
                    print(x.getDjangoModel())

                for x in getSubclassesFromModule(bases.DescriptorCalculator, descriptors):
                    print(x.getDjangoModel())
