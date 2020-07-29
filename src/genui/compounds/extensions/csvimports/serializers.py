"""
serializers

Created by: Martin Sicho
On: 7/16/20, 1:52 PM
"""
from . import models
from genui.compounds.extensions.fileimports.serializer import FileSetSerializer, FileSetUpdateSerializer

class CSVSetSerializer(FileSetSerializer):

    class Meta:
        model = models.CSVCompounds
        file_extension = '.csv'
        fields = FileSetSerializer.Meta.fields + ('nameCol', 'smilesCol', 'activityCol', 'activityTypeCol', 'activityUnitsCol', 'colSeparator', 'emptyValue')
        read_only_fields = FileSetSerializer.Meta.read_only_fields

class CSVSetUpdateSerializer(FileSetUpdateSerializer):

    class Meta:
        model = models.CSVCompounds
        fields = FileSetUpdateSerializer.Meta.fields
        read_only_fields = FileSetUpdateSerializer.Meta.read_only_fields
