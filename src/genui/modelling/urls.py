"""
urls

Created by: Martin Sicho
On: 5/3/20, 6:13 PM
"""
from genui.modelling.apps import ModellingConfig
from genui.utils.inspection import discover_extensions_urlpatterns

urlpatterns = discover_extensions_urlpatterns(ModellingConfig.name)

