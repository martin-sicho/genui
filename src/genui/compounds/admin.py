from django.contrib import admin

# Register your models here.
from genui.compounds.models import MolSet, ActivitySet, Molecule

@admin.register(Molecule)
class MolAdmin(admin.ModelAdmin):
    pass

@admin.register(MolSet)
class MolSetAdmin(admin.ModelAdmin):
    readonly_fields = ["created", "updated"]

@admin.register(ActivitySet)
class ActivitySetAdmin(admin.ModelAdmin):
    readonly_fields = ["created", "updated"]

