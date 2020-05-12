from django.contrib import admin
from django.contrib.auth.models import Permission


@admin.register(Permission)
class PermissionsAdmin(admin.ModelAdmin):
    pass
