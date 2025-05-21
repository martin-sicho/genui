from rest_framework import serializers
from . import models
from genui.compounds.serializers import MolSetSerializer, MolSetUpdateSerializer

class PapyrusAssaySerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.PapyrusAssay
        fields = ('assayID',)

class PapyrusTargetSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.PapyrusTarget
        fields = ('targetID',)


class PapyrusActivitySerializer(serializers.ModelSerializer):

    class Meta:
        model = models.PapyrusActivity
        fields = (
            'value',
            'type',
            'relation',
            'assay',
            'target',
        )

class PapyrusSetSerializer(MolSetSerializer):
    targets = PapyrusTargetSerializer(many=True)
    activities = serializers.PrimaryKeyRelatedField(many=True, read_only=True)

    class Meta:
        model = models.PapyrusCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'targets', 'activities')
        read_only_fields = ('created', 'updated')

class PapyrusMoleculeSerializer(serializers.ModelSerializer):

    class Meta:
        model = models.PapyrusMolecule
        fields = '__all__'

class PapyrusSetInitSerializer(PapyrusSetSerializer):
    maxPerTarget = serializers.IntegerField(min_value=1, required=False)
    taskID = serializers.CharField(required=False, read_only=True)
    targets = serializers.ListField(child=serializers.CharField(), min_length=1, required=True, write_only=True)

    def create(self, validated_data):
        instance = models.PapyrusCompounds(
            name=validated_data["name"]
            , description=validated_data["description"]
            , project=validated_data["project"]
        )
        instance.save()
        for target in validated_data['targets']:
            target = models.PapyrusTarget.objects.get_or_create(targetID=target)[0]
            instance.targets.add(target)
        instance.save()
        return instance

    class Meta:
        model = models.PapyrusCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'maxPerTarget', 'taskID', 'targets', 'activities')
        read_only_fields = ('created', 'updated', 'taskID', 'targets', 'activities')

class PapyrusSetUpdateSerializer(MolSetUpdateSerializer):

    class Meta:
        model = models.PapyrusCompounds
        fields = MolSetUpdateSerializer.Meta.fields
        read_only_fields = MolSetUpdateSerializer.Meta.read_only_fields
