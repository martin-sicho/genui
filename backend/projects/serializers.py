"""
serializers

Created by: Martin Sicho
On: 05-12-19, 12:25
"""
from django.conf import settings

from modelling import helpers
from .models import Project
from rest_framework import serializers

# Serializers define the API representation.
class ProjectSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Project
        fields = ('id', 'name', 'description', 'created', 'updated')
        read_only_fields = ('created', 'updated')

    def create(self, validated_data):
        ret = super().create(validated_data)

        for app in settings.GENUI_MODEL_APPS:
            helpers.createDefaultModels(ret, app)

        return ret

