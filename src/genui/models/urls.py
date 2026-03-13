from genui.models.apps import ModelsConfig
from genui.utils.inspection import discover_extensions_urlpatterns

urlpatterns = discover_extensions_urlpatterns(ModelsConfig.name)

