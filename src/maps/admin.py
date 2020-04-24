from django.contrib import admin
from . import models

@admin.register(models.Map)
class MapAdmin(admin.ModelAdmin):
    pass
