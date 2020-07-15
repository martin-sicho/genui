"""
serializers

Created by: Martin Sicho
On: 7/13/20, 3:58 PM
"""
from rest_framework import serializers

from genui.compounds.serializers import MolSetSerializer
from . import models
from genui.compounds.models import MolSetFile
import uuid
import re


class SDFSetSerializer(MolSetSerializer):
    file = serializers.FileField(use_url=True)

    class Meta:
        model = models.SDFCompounds
        fields = MolSetSerializer.Meta.fields + ('activitiesProp', 'activityTypesProp', 'activityUnitsProp', 'dataSeparator', 'file')
        read_only_fields = MolSetSerializer.Meta.read_only_fields

    def create(self, validated_data):
        uploaded_file = validated_data['file']
        del validated_data['file']
        instance =  super().create(validated_data)
        MolSetFile.create(molset=instance, filename=f'{re.sub("[^0-9a-zA-Z]+", "_", instance.name)}_{uuid.uuid4().hex}.sdf', file=uploaded_file)

        return instance

