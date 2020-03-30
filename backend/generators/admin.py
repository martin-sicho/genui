from django.contrib import admin
from . import models

@admin.register(models.Generator)
class GeneratorAdmin(admin.ModelAdmin):
    pass

@admin.register(models.DrugExNet)
class DrugExNetAdmin(admin.ModelAdmin):
    pass

@admin.register(models.DrugExAgent)
class DrugExAgentAdmin(admin.ModelAdmin):
    pass
