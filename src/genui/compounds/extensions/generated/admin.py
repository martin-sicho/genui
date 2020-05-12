from django.contrib import admin

from . import models
from genui.compounds.admin import MolSetAdmin


@admin.register(models.GeneratedMolSet)
class ChEMBLCompoundsAdmin(MolSetAdmin):
    pass
