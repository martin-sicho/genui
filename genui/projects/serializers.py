"""
serializers

Created by: Martin Sicho
On: 05-12-19, 12:25
"""

from .models import Project
from rest_framework import serializers

# Serializers define the API representation.
class GenUIProjectSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Project
        fields = ('id', 'name', 'description', 'created', 'updated')
        read_only_fields = ('created', 'updated')

