from django.contrib import admin

from genui.compounds.admin import MolSetAdmin, ActivitySetAdmin
from . import models


@admin.register(models.ChEMBLCompounds)
class ChEMBLCompoundsAdmin(MolSetAdmin):
    pass

@admin.register(models.ChEMBLActivities)
class ChEMBLActivitiesAdmin(ActivitySetAdmin):
    pass
