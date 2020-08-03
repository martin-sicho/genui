"""
serializers

Created by: Martin Sicho
On: 25-02-20, 16:35
"""
from rest_framework import serializers

from genui.compounds.models import MolSet, Molecule
from genui.compounds.serializers import GenericMolSetSerializer, MoleculeSerializer
from genui.models.serializers import ModelSerializer, TrainingStrategySerializer, TrainingStrategyInitSerializer
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
        fields = TrainingStrategySerializer.Meta.fields + ("descriptors",)

class MapSerializer(ModelSerializer):
    trainingStrategy = MappingStrategySerializer(many=False)
    molsets = GenericMolSetSerializer(many=True, required=True, allow_null=False)

    class Meta:
        model = models.Map
        fields = [x for x in ModelSerializer.Meta.fields if x not in ('validationStrategy', 'performance')] + ['molsets']
        read_only_fields = [x for x in ModelSerializer.Meta.read_only_fields if x not in ('validationStrategy', 'performance')]

class MapInitSerializer(MapSerializer):
    trainingStrategy = MappingStrategyInitSerializer(many=False)
    molsets = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all(), required=True, allow_null=False, allow_empty=False)

    class Meta:
        model = models.Map
        fields = MapSerializer.Meta.fields
        read_only_fields = MapSerializer.Meta.read_only_fields

    def create(self, validated_data, **kwargs):
        instance = super().create(validated_data, **kwargs)
        instance.molsets.set(validated_data['molsets'])
        instance.save()

        ts_data = validated_data['trainingStrategy']
        ts = models.MappingStrategy.objects.create(
            modelInstance=instance,
            algorithm = ts_data['algorithm'],
            mode = ts_data['mode'],
        )
        ts.descriptors.set(ts_data['descriptors'])
        ts.save()
        self.saveParameters(ts, ts_data)

        return instance


class PointSerializer(serializers.ModelSerializer):
    molecule = serializers.PrimaryKeyRelatedField(many=False, queryset=Molecule.objects.all())
    compoundSets = serializers.PrimaryKeyRelatedField(many=True, queryset=MolSet.objects.all())
    map = serializers.PrimaryKeyRelatedField(many=False, queryset=models.Map.objects.all())

    class Meta:
        model = models.Point
        fields = ('id', 'x', 'y', 'map', 'molecule', 'compoundSets')