from django.contrib import admin

from .models import GenUIProject

@admin.register(GenUIProject)
class GenUIProjectAdmin(admin.ModelAdmin):

    readonly_fields = ["created", "updated"]
