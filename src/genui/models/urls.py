"""
urls

Created by: Martin Sicho
On: 5/3/20, 6:13 PM
"""
from genui.models.apps import ModelsConfig
from genui.utils.inspection import discover_extensions_urlpatterns

urlpatterns = discover_extensions_urlpatterns(ModelsConfig.name)

