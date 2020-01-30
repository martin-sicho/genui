from django.apps import AppConfig


class QsarConfig(AppConfig):
    name = 'qsar'

    def ready(self, force_inspect=False):
        from modelling import helpers
        from .core import bases
        helpers.inspectCore("qsar", force=force_inspect, modules=["builders", "descriptors"], additional_bases=[bases.DescriptorCalculator])
