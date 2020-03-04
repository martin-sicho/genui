from django.contrib import admin
from . import models

@admin.register(models.QSARModel)
class QSARModelAdmin(admin.ModelAdmin):
    readonly_fields = ["created", "updated"]
