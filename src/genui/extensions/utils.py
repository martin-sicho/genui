"""
utils

Created by: Martin Sicho
On: 4/27/20, 6:35 PM
"""
import importlib
import sys

from genui import apps
from genui.commons.inspection import importFromPackage

from django.urls import path, include

def discover_extensions(packages):
    ret = []
    for package_name in packages:
        try:
            package = importlib.import_module(package_name)
        except ModuleNotFoundError as exp:
            sys.stderr.write(f"WARNING: Failed to find extensions module: {package_name}. Skipping...\n")
            continue

        subs = package.__all__
        for sub in subs:
            ret.append(f"{package_name}.{sub}")
    return ret

def disover_app_urls_module(app_name, parent=None):
    app = importlib.import_module(app_name)

    error_skipped = False
    if not importFromPackage(app, 'genuisetup', exception=False):
        sys.stderr.write(f"WARNING: {app_name}.genuisetup module not found.\n")
        error_skipped = True
    elif not importFromPackage(app, 'urls', exception=False):
        sys.stderr.write(f"WARNING: No {app_name}.urls module found for app:  {app_name}\n")
        error_skipped = True
    elif not hasattr(app.genuisetup, 'PARENT'):
        if parent == 'genui':
            return app.urls
        else:
            return None
    elif app.genuisetup.PARENT != parent:
        return None

    if error_skipped:
        sys.stderr.write(f"WARNING: Skipping urls discovery for: {app_name}\n")
        return None

    return app.urls

def discover_apps_urls(app_names, prefix='', app_names_as_root=False):
    urls = []
    for app in app_names:
        urls_module = disover_app_urls_module(app, parent='genui')
        if urls_module:
            urls.append(path(
                    f'{prefix.rstrip("/") + "/" if prefix else ""}{app.split(".")[1] + "/" if app_names_as_root else ""}',
                    include(urls_module.__name__)
                )
            )

    return urls

def discover_extensions_urlpatterns(parent):
    ret = []
    for extension_name in apps.extensions():
        urls = disover_app_urls_module(extension_name, parent)
        if urls:
            ret.extend(urls.urlpatterns)
    return ret

