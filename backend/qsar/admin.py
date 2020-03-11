from django.contrib import admin
from . import models

class UnitsInline(admin.TabularInline):
    model = models.QSARTrainingStrategy

@admin.register(models.QSARModel)
class QSARModelAdmin(admin.ModelAdmin):
    readonly_fields = ["created", "updated"]
    inlines = [
        UnitsInline,
    ]
