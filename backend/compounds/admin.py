from django.contrib import admin

# Register your models here.
from compounds.models import ChEMBLCompounds, ChEMBLActivities, MolSet, ActivitySet


@admin.register(MolSet)
class MolSetAdmin(admin.ModelAdmin):
    readonly_fields = ["created", "updated"]

@admin.register(ChEMBLCompounds)
class ChEMBLCompoundsAdmin(MolSetAdmin):
    pass

@admin.register(ActivitySet)
class ActivitySetAdmin(admin.ModelAdmin):
    readonly_fields = ["created", "updated"]

@admin.register(ChEMBLActivities)
class ChEMBLActivitiesAdmin(ActivitySetAdmin):
    pass

