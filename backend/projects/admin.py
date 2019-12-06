from django.contrib import admin

from .models import Project

@admin.register(Project)
class GenUIProjectAdmin(admin.ModelAdmin):

    readonly_fields = ["created", "updated"]
