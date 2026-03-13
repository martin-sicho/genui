from genui.utils.inspection import discover_extensions_urlpatterns
from genui.utils.apps import UtilsConfig

urlpatterns = discover_extensions_urlpatterns(UtilsConfig.name)
