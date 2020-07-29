"""
serializers

Created by: Martin Sicho
On: 4/27/20, 8:47 PM
"""
from rest_framework import serializers
from . import models
from genui.compounds.serializers import MolSetSerializer, MolSetUpdateSerializer

class ChEMBLAssaySerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.ChEMBLAssay
        fields = ('assayID',)

class ChEMBLTargetSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.ChEMBLTarget
        fields = ('targetID',)

class ChEMBLSetSerializer(MolSetSerializer):
    targets = ChEMBLTargetSerializer(many=True)

    class Meta:
        model = models.ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'targets', 'activities')
        read_only_fields = ('created', 'updated')

class ChEMBLSetInitSerializer(ChEMBLSetSerializer):
    maxPerTarget = serializers.IntegerField(min_value=1, required=False)
    taskID = serializers.CharField(required=False, read_only=True)
    targets = serializers.ListField(child=serializers.CharField(), min_length=1, required=True, write_only=True)

    def create(self, validated_data):
        instance = models.ChEMBLCompounds(
            name=validated_data["name"]
            , description=validated_data["description"]
            , project=validated_data["project"]
        )
        instance.save()
        for target in validated_data['targets']:
            target = models.ChEMBLTarget.objects.get_or_create(targetID=target)[0]
            instance.targets.add(target)
        instance.save()
        return instance

    class Meta:
        model = models.ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'maxPerTarget', 'taskID', 'targets', 'activities')
        read_only_fields = ('created', 'updated', 'taskID', 'targets', 'activities')

class ChEMBLSetUpdateSerializer(MolSetUpdateSerializer):

    class Meta:
        model = models.ChEMBLCompounds
        fields = MolSetUpdateSerializer.Meta.fields
        read_only_fields = MolSetUpdateSerializer.Meta.read_only_fields
