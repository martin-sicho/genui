import sys

from django.apps import AppConfig

class QsarConfig(AppConfig):
    name = 'qsar'

    def ready(self):
        if sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate'):
            from .algorithms import bases
            from commons.helpers import getSubclassesFromModule
            from .algorithms import algorithms as algs
            from .algorithms import metrics as metrics

            for x in getSubclassesFromModule(bases.BaseAlgorithm, algs):
                x.getParams()
                x.getDjangoModel()

            for x in getSubclassesFromModule(bases.ValidationMetric, metrics):
                x.getDjangoModel()
