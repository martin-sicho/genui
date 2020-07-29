"""
serializers

Created by: Martin Sicho
On: 7/13/20, 3:58 PM
"""
from . import models
from genui.compounds.extensions.fileimports.serializer import FileSetSerializer, FileSetUpdateSerializer


class SDFSetSerializer(FileSetSerializer):

    class Meta:
        model = models.SDFCompounds
        file_extension = '.sdf'
        fields = FileSetSerializer.Meta.fields + ('activitiesProp', 'activityTypesProp', 'activityUnitsProp', 'dataSeparator')
        read_only_fields = FileSetSerializer.Meta.read_only_fields

class SDFSetUpdateSerializer(FileSetUpdateSerializer):

    class Meta:
        model = models.SDFCompounds
        fields = FileSetUpdateSerializer.Meta.fields
        read_only_fields = FileSetUpdateSerializer.Meta.read_only_fields

