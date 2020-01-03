"""
serializers

Created by: Martin Sicho
On: 18-12-19, 10:27
"""
from rest_framework import serializers

from projects.models import Project
from .models import MolSet, Molecule, ChEMBLCompounds, ChEMBLTarget, ChEMBLAssay, ChEMBLActivities


class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    providers = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all())

    class Meta:
        model = Molecule
        fields = ('id', 'canonicalSMILES', 'inchiKey', 'providers')

class ChEMBLAssaySerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ChEMBLAssay
        fields = ('assayID',)

class ChEMBLTargetSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ChEMBLTarget
        fields = ('targetID',)

class MolSetSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())

    class AutoSchemaMixIn:
        def get_operation(self, path, method):
            ret = super().get_operation(path, method)
            if method in ('POST', 'PUT', 'PATCH'):
                ret['responses']['200']['content']['application/json']['schema']['properties']['task'] = {
                    'type' : 'string'
                }
            return ret

    class Meta:
        model = MolSet
        fields = ('id', 'name', 'description', 'created', 'updated', 'project')
        read_only_fields = ('created', 'updated')

class GenericMolSetSerializer(MolSetSerializer):
    className = serializers.CharField(default="")

    class Meta:
        model = MolSet
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'className')
        read_only_fields = ('created', 'updated')

    def to_representation(self, instance):
        ret = super().to_representation(instance)
        ret['className'] = instance.__class__.__name__
        return ret

class ChEMBLSetSerializer(MolSetSerializer):
    targets = ChEMBLTargetSerializer(many=True)
    activities = serializers.PrimaryKeyRelatedField(many=False, read_only=True) #TODO: write a proper serializer for the activities

    class Meta:
        model = ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'targets', 'activities')
        read_only_fields = ('created', 'updated', 'activities')

    def update(self, instance, validated_data):
        super().update(instance, validated_data)
        # FIXME: this needs to be implemented in order for PUT and PATCH to work

class ChEMBLSetInitSerializer(ChEMBLSetSerializer):
    maxPerTarget = serializers.IntegerField(min_value=1, required=False)
    taskID = serializers.CharField(required=False, read_only=True)
    targets = serializers.ListField(child=serializers.CharField(), min_length=1, required=True, write_only=True)

    def create(self, validated_data):
        instance = ChEMBLCompounds(
            name=validated_data["name"]
            , description=validated_data["description"]
            , project=validated_data["project"]
        )
        instance.save()
        for target in validated_data['targets']:
            target = ChEMBLTarget.objects.get_or_create(targetID=target)[0]
            instance.targets.add(target)
        instance.save()
        return instance

    class Meta:
        model = ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'targets', 'maxPerTarget', 'taskID', 'targets', 'activities')
        read_only_fields = ('created', 'updated', 'taskID', 'targets', 'activities')