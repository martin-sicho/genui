"""
serializers

Created by: Martin Sicho
On: 7/13/20, 3:58 PM
"""
from . import models
from genui.compounds.extensions.fileimports.serializer import FileSetSerializer


class SDFSetSerializer(FileSetSerializer):

    class Meta:
        model = models.SDFCompounds
        file_extension = '.sdf'
        fields = FileSetSerializer.Meta.fields + ('activitiesProp', 'activityTypesProp', 'activityUnitsProp', 'dataSeparator')
        read_only_fields = FileSetSerializer.Meta.read_only_fields

