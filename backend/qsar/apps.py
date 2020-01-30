from django.apps import AppConfig


class QsarConfig(AppConfig):
    name = 'qsar'

    def ready(self, force_inspect=False):
        from modelling import helpers
        helpers.inspectCore("qsar", force=force_inspect)
