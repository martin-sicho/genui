"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 9:34 AM
"""
import importlib

from django.core.management import BaseCommand


class Command(BaseCommand):
    help = 'Sets up the genui app and the extensions.'

    def add_arguments(self, parser):
        parser.add_argument(
            '--force',
            action='store_true',
            help='Force all updates',
        )

        parser.add_argument(
            '--strict',
            action='store_true',
            help='Do not ignore some exceptions. Such as when the genuisetup module is not found.',
        )

    def handle(self, *args, **options):
        apps = []
        try:
            from django.conf import settings
            apps = settings.GENUI_SETTINGS['APPS']
        except Exception as exp:
            # FIXME: make this not catch-all
            self.style.WARNING('Failed to load GENUI_SETTINGS from settings.py. Loading internal modules only...')
            from genui import apps
            apps = apps.all_()

        for app in apps:
            try:
                setupmodule = importlib.import_module(f'{app}.genuisetup')
                setupmodule.setup(
                    force=bool(options['force']),
                    strict=bool(options['strict'])
                )
            except ModuleNotFoundError as exp:
                if not options['strict']:
                    self.stderr.write(self.style.WARNING(f'Failed to find the genuisetup module for app or extension: "{app}". No setup will be done. Reason: {exp}'))
                    continue
                else:
                    raise exp
            self.stdout.write(self.style.SUCCESS('Successful setup for: "%s"' % app))