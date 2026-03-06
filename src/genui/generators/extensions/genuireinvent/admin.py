from django.contrib import admin
from . import models

@admin.register(models.ReinventNet)
class ReinventNetAdmin(admin.ModelAdmin):
    pass

@admin.register(models.ReinventAgent)
class ReinventAgentAdmin(admin.ModelAdmin):
    pass