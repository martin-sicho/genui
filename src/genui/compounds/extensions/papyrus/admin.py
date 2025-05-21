from django.contrib import admin

from genui.compounds.admin import MolSetAdmin, ActivitySetAdmin
from . import models


@admin.register(models.PapyrusCompounds)
class PapyrusCompoundsAdmin(MolSetAdmin):
    pass

@admin.register(models.PapyrusActivities)
class PapyrusActivitiesAdmin(ActivitySetAdmin):
    pass
