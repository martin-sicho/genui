"""
urls

Created by: Martin Sicho
On: 5/3/20, 5:34 PM
"""
from genui.extensions.utils import discover_extensions_urlpatterns
from genui.utils.apps import UtilsConfig

urlpatterns = discover_extensions_urlpatterns(UtilsConfig.name)