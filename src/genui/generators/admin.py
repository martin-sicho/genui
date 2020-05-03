from django.contrib import admin

import genui.generators.extensions.genuidrugex.models
from . import models

@admin.register(models.Generator)
class GeneratorAdmin(admin.ModelAdmin):
    pass

@admin.register(genui.generators.extensions.genuidrugex.models.DrugExNet)
class DrugExNetAdmin(admin.ModelAdmin):
    pass

@admin.register(genui.generators.extensions.genuidrugex.models.DrugExAgent)
class DrugExAgentAdmin(admin.ModelAdmin):
    pass
