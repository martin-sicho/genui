from django.contrib import admin
from . import models

@admin.register(models.Generator)
class GeneratorAdmin(admin.ModelAdmin):
    pass
