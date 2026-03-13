from rest_framework import serializers

from genui.compounds.models import MolSet, Molecule
from genui.compounds.serializers import GenericMolSetSerializer
from genui.models.serializers import ModelSerializer, TrainingStrategySerializer, TrainingStrategyInitSerializer, \
    ModelFileSerializer
from genui.qsar.models import EmbeddingCalculator
from genui.qsar.serializers import EmbeddingCalculatorSerializer
from . import models

class MappingStrategySerializer(TrainingStrategySerializer):
    embeddings = EmbeddingCalculatorSerializer(many=True)

    class Meta:
        model = models.MappingStrategy
        fields = TrainingStrategySerializer.Meta.fields + ("embeddings",)

class MappingStrategyInitSerializer(TrainingStrategyInitSerializer):
    embeddings = serializers.PrimaryKeyRelatedField(many=True, queryset=EmbeddingCalculator.objects.all(), allow_empty=False)

    class Meta:
        model = models.MappingStrategy
        fields = TrainingStrategyInitSerializer.Meta.fields + ("embeddings",)

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

    def is_valid(self, *, raise_exception=False):
        initial_data = self.initial_data
        if "trainingStrategy" in initial_data and "embeddings" in initial_data["trainingStrategy"]:
            embeddings = []
            for emb in initial_data["trainingStrategy"]["embeddings"]:
                if isinstance(emb, dict):
                    eid, _ = EmbeddingCalculator.objects.get_or_create(**emb)
                    embeddings.append(eid.id)
                else:
                    embeddings.append(emb)
            initial_data["trainingStrategy"]["embeddings"] = embeddings
        return super().is_valid(raise_exception=raise_exception)

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
        ts.embeddings.set(ts_data['embeddings'])
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
