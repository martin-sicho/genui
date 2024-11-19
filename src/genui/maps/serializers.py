"""
serializers

Created by: Martin Sicho
On: 25-02-20, 16:35
"""
from rest_framework import serializers

from genui.compounds.models import MolSet, Molecule
from genui.compounds.serializers import GenericMolSetSerializer, MoleculeSerializer
from genui.models.serializers import ModelSerializer, TrainingStrategySerializer, TrainingStrategyInitSerializer, \
    ModelFileSerializer
from genui.qsar.serializers import DescriptorGroupSerializer
from . import models

class MappingStrategySerializer(TrainingStrategySerializer):
    descriptors = DescriptorGroupSerializer(many=True)

    class Meta:
        model = models.MappingStrategy
        fields = TrainingStrategySerializer.Meta.fields + ("descriptors",)

class MappingStrategyInitSerializer(TrainingStrategyInitSerializer):
    descriptors = serializers.PrimaryKeyRelatedField(many=True, queryset=models.DescriptorGroup.objects.all(), allow_empty=False)

    class Meta:
        model = models.MappingStrategy
        fields = TrainingStrategyInitSerializer.Meta.fields + ("descriptors",)

class MapSerializer(ModelSerializer):
    trainingStrategy = MappingStrategySerializer(many=False)
    molsets = GenericMolSetSerializer(many=True, required=True, allow_null=False)
    chemspaceJSON = ModelFileSerializer(allow_null=True, read_only=True)

    class Meta:
        model = models.Map
        fields = ModelSerializer.Meta.fields + ('molsets', 'chemspaceJSON')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('chemspaceJSON',)

class MapInitSerializer(MapSerializer):
    trainingStrategy = MappingStrategyInitSerializer(many=False)
    molsets = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all(), required=True, allow_null=False, allow_empty=False)

    class Meta:
        model = models.Map
        fields = MapSerializer.Meta.fields
        read_only_fields = MapSerializer.Meta.read_only_fields

    def create(self, validated_data, **kwargs):
        molsets = validated_data.pop('molsets')
        ts_data = validated_data.pop('trainingStrategy')
        
        # Create the instance using the parent's create method
        instance = super().create(validated_data, **kwargs)
        
        # Set the molsets
        instance.molsets.set(molsets)
        
        # Create and set the training strategy
        ts = models.MappingStrategy.objects.create(
            modelInstance=instance,
            algorithm=ts_data['algorithm'],
            mode=ts_data['mode'],
        )
        ts.descriptors.set(ts_data['descriptors'])
        ts.save()
        
        # Save parameters
        self.saveParameters(ts, ts_data)
        
        return instance


class PointSerializer(serializers.ModelSerializer):
    molecule = serializers.PrimaryKeyRelatedField(many=False, queryset=Molecule.objects.all())
    compoundSets = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all())
    map = serializers.PrimaryKeyRelatedField(many=False, queryset=models.Map.objects.all())

    class Meta:
        model = models.Point
        fields = ('id', 'x', 'y', 'map', 'molecule', 'smiles', 'compoundSets')