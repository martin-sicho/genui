"""
serializers

Created by: Martin Sicho
On: 18-12-19, 10:27
"""
from rest_framework import serializers
from rest_framework.schemas.openapi import AutoSchema

from projects.models import Project
from .models import MolSet, Molecule, ChEMBLCompounds, ChEMBLTarget, ChEMBLAssay

class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    providers = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all())

    class Meta:
        model = Molecule
        fields = ('id', 'canonicalSMILES', 'inchiKey', 'providers')

class ChEMBLAssaySerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ChEMBLAssay
        fields = ('id', 'assayID',)
        read_only_fields = ('id',)
        extra_kwargs = {
            'assayID': {
                'validators': [],
            },
        }

class ChEMBLTargetSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ChEMBLTarget
        fields = ('id', 'targetID',)
        read_only_fields = ('id',)
        extra_kwargs = {
            'targetID': {
                'validators': [],
            },
        }

class MolSetSerializer(serializers.HyperlinkedModelSerializer):
    project = serializers.PrimaryKeyRelatedField(many=False, queryset=Project.objects.all())

    class Schema(AutoSchema):
        def get_operation(self, path, method):
            ret = super().get_operation(path, method)
            if method in ('POST', 'PUT', 'PATCH'):
                ret['responses']['200']['content']['application/json']['schema']['properties']['task'] = {
                    'type' : 'string'
                }
            return ret

    class Meta:
        model = MolSet

class ChEMBLSetSerializer(MolSetSerializer):
    assays = ChEMBLAssaySerializer(many=True, required=False)
    targets = ChEMBLTargetSerializer(many=True, required=False)
    # TODO: add validation that checks at least one of the fields is supplied

    class Meta:
        model = ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'assays', 'targets')
        read_only_fields = ('created', 'updated')

    def create(self, validated_data):
        targets = []
        if 'targets' in validated_data:
            for target in validated_data['targets']:
                targets.append(ChEMBLTarget.objects.get_or_create(targetID=target['targetID'])[0])

        assays = []
        if 'assays' in validated_data:
            for assay in validated_data['assays']:
                assays.append(ChEMBLAssay.objects.get_or_create(assayID=assay['assayID'])[0])

        instance = ChEMBLCompounds(
            name=validated_data["name"]
            , description=validated_data["description"]
            , project=validated_data["project"]
        )
        instance.save()
        instance.assays.add(*assays)
        instance.targets.add(*targets)
        instance.save()
        return instance

    def update(self, instance, validated_data):
        super().update(instance, validated_data)
        # FIXME: this needs to be implemented in order for PUT and PATCH to work

class TaskSerializer(serializers.Serializer):
    task_id = serializers.CharField(allow_blank=False)
    status = serializers.CharField(allow_blank=False)

class MolSetTasksSerializerFactory:
    class Schema(AutoSchema):
        def get_operation(self, path, method):
            ret = super().get_operation(path, method)
            ret['responses']['200']['content']['application/json']['schema'] = {
                'type' : 'object',
                'properties' : {
                        'taskName' : {
                            'type' : 'array',
                            'items' : {
                                'properties' : {
                                    'task_id' : {
                                        'type' : 'string'
                                    },
                                    'status' : {
                                        'type' : 'string'
                                    },
                                    'required': ['task_id', 'status']
                                }
                            }
                        },
                        'required': ['taskName']
                    }
                }
            return ret

    @staticmethod
    def get(field_names):
        return type('MolSetTasksSerializer', (serializers.Serializer,), {x : serializers.ListField(child=TaskSerializer(required=False), allow_empty=True) for x in field_names})