from django.contrib import admin

# Register your models here.
from compounds.models import ChEMBLCompounds


@admin.register(ChEMBLCompounds)
class ChEMBLCompoundsAdmin(admin.ModelAdmin):

    readonly_fields = ["created", "updated"]
