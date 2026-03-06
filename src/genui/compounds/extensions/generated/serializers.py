"""
serializers

Created by: Martin Sicho
On: 5/12/20, 9:42 AM
"""
from rest_framework import serializers

from genui.compounds.serializers import GenericMolSetSerializer, MolSetUpdateSerializer
from genui.generators.serializers import GeneratorSerializer
from genui.projects.models import Project
from . import models


class GeneratedSetSerializer(GenericMolSetSerializer):
    source = GeneratorSerializer(many=False)

    class Meta:
        model = models.GeneratedMolSet
        fields = list(GenericMolSetSerializer.Meta.fields) + ['source']
        read_only_fields = list(GenericMolSetSerializer.Meta.read_only_fields)


class GeneratedSetInitSerializer(GeneratedSetSerializer):
    source = serializers.PrimaryKeyRelatedField(many=False, queryset=models.Generator.objects.all())
    nSamples = serializers.IntegerField(min_value=1, default=1000)
    minScore = serializers.FloatField(required=False, allow_null=True, default=None)
    taskID = serializers.UUIDField(required=False, read_only=True)

    class Meta:
        model = models.GeneratedMolSet
        fields = GeneratedSetSerializer.Meta.fields + ["nSamples", "minScore", "taskID"]
        read_only_fields = GeneratedSetSerializer.Meta.read_only_fields + ["taskID"]

    def create(self, validated_data):
        molset_class = models.GeneratedMolSet
        if "className" in validated_data:
            # TODO:  make this more general (look in different modules as well)
            molset_class = getattr(models, validated_data["className"])

        # Get extraArgs (but exclude minScore - it's not a DB field)
        extraArgs = validated_data.get("extraArgs", {}).copy()

        # Extract minScore separately (don't include in extraArgs)
        min_score = validated_data.get("minScore")

        # Create instance (extraArgs should only contain valid DB fields)
        instance = molset_class.objects.create(
            name=validated_data["name"],
            description=validated_data["description"] if "description" in validated_data else "",
            source=validated_data["source"],
            project=validated_data["project"],
            **extraArgs
        )

        # Set nSamples as instance attribute (used by task)
        instance.nSamples = validated_data["nSamples"]

        # Set minScore as instance attribute (used by initializer, not saved to DB)
        if min_score is not None:
            instance.minScore = min_score

        return instance

class GeneratedSetUpdateSerializer(MolSetUpdateSerializer):

    class Meta:
        model = models.GeneratedMolSet
        fields = MolSetUpdateSerializer.Meta.fields
        read_only_fields = MolSetUpdateSerializer.Meta.read_only_fields