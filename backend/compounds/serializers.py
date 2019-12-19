"""
serializers

Created by: Martin Sicho
On: 18-12-19, 10:27
"""
from rest_framework import serializers

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
    populateTask = serializers.CharField(read_only=True, allow_null=True)

    class Meta:
        model = MolSet

class ChEMBLSetSerializer(MolSetSerializer):
    assays = ChEMBLAssaySerializer(many=True, required=False)
    targets = ChEMBLTargetSerializer(many=True, required=False)
    # TODO: add validation that checks at least one of the fields is supplied

    class Meta:
        model = ChEMBLCompounds
        fields = ('id', 'name', 'description', 'created', 'updated', 'project', 'assays', 'targets', 'populateTask')
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
        # FIXME: this should probably be implemented

