"""
serializer

Created by: Martin Sicho
On: 7/16/20, 2:04 PM
"""
from rest_framework import serializers
import uuid
import re

from genui.compounds.models import MolSetFile
from genui.compounds.serializers import MolSetSerializer


class FileSetSerializer(MolSetSerializer):
    file = serializers.FileField(use_url=True)

    class Meta:
        fields = MolSetSerializer.Meta.fields + ('file',)
        read_only_fields = MolSetSerializer.Meta.read_only_fields
        file_extension = ''

    def create(self, validated_data):
        uploaded_file = validated_data['file']
        del validated_data['file']
        instance =  super().create(validated_data)
        MolSetFile.create(molset=instance, filename=f'{re.sub("[^0-9a-zA-Z]+", "_", instance.name)}_{uuid.uuid4().hex}{self.Meta.file_extension}', file=uploaded_file)

        return instance

