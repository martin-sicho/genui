from django.apps import AppConfig
import sys
from django.db import transaction

class ModellingConfig(AppConfig):
    name = 'modelling'

    def ready(self):
        from . import signals
        if len(sys.argv) > 1 and sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate'):
            from .core import bases
            from commons.helpers import getSubclassesFromModule
            from .core import algorithms as algs
            from .core import metrics

            with transaction.atomic():
                for x in getSubclassesFromModule(bases.Algorithm, algs):
                    print(x.getParams())
                    print(x.getDjangoModel())

                for x in getSubclassesFromModule(bases.ValidationMetric, metrics):
                    print(x.getDjangoModel())
