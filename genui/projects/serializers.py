"""
serializers

Created by: Martin Sicho
On: 05-12-19, 12:25
"""

from .models import GenUIProject
from rest_framework import serializers

# Serializers define the API representation.
class GenUIProjectSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = GenUIProject
        fields = ('id', 'name', 'description', 'created', 'updated')
        read_only_fields = ('created', 'updated')

