from django.contrib import admin

# Register your models here.
from compounds.models import ChEMBLCompounds, ChEMBLActivities


@admin.register(ChEMBLCompounds)
class ChEMBLCompoundsAdmin(admin.ModelAdmin):

    readonly_fields = ["created", "updated"]

@admin.register(ChEMBLActivities)
class ChEMBLActivitiesAdmin(admin.ModelAdmin):

    readonly_fields = ["created", "updated"]
