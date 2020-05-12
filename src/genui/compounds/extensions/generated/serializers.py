"""
serializers

Created by: Martin Sicho
On: 5/12/20, 9:42 AM
"""
from rest_framework import serializers

from genui.compounds.serializers import GenericMolSetSerializer
from genui.generators.serializers import GeneratorSerializer
from . import models


class GeneratedSetSerializer(GenericMolSetSerializer):
    source = GeneratorSerializer(many=False)

    class Meta:
        model = models.GeneratedMolSet
        fields = list(GenericMolSetSerializer.Meta.fields) + ['source']
        read_only_fields = list(GenericMolSetSerializer.Meta.read_only_fields)


class GeneratedSetInitSerializer(GeneratedSetSerializer):
    source = serializers.PrimaryKeyRelatedField(many=False, queryset=models.Generator.objects.all())
    nSamples = serializers.IntegerField(min_value=1)
    taskID = serializers.UUIDField(required=False, read_only=True)

    class Meta:
        model = models.GeneratedMolSet
        fields = GeneratedSetSerializer.Meta.fields + ["nSamples", "taskID"]
        read_only_fields = GeneratedSetSerializer.Meta.read_only_fields + ["taskID"]

    def create(self, validated_data):
        molset_class = models.GeneratedMolSet
        if "className" in validated_data:
            # TODO:  make this more general (look in different modules as well)
            molset_class = getattr(models, validated_data["className"])
        extraArgs = validated_data["extraArgs"] if "extraArgs" in validated_data else dict()
        instance = molset_class.objects.create(
            name=validated_data["name"],
            description=validated_data["description"] if "description" in validated_data else "",
            source=validated_data["source"],
            project=validated_data["project"],
            **extraArgs

        )
        instance.nSamples = validated_data["nSamples"]

        return instance