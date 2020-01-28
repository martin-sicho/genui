"""
serializers

Created by: Martin Sicho
On: 27-01-20, 17:00
"""
from rest_framework import serializers

from compounds.models import MolSet
from compounds.serializers import MolSetSerializer
from modelling.models import ModelPerformanceMetric
from modelling.serializers import ModelSerializer, ValidationStrategySerializer, TrainingStrategySerializer, \
    TrainingStrategyInitSerializer
from . import models

class DrugExCorpusSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DrugeExCorpus
        fields = ("corpusFile", "vocFile")
        read_only_fields = ("corpusFile", "vocFile")

class DrugExValidationStrategySerializer(ValidationStrategySerializer):

    class Meta:
        model = models.DrugExValidationStrategy
        fields = ValidationStrategySerializer.Meta.fields + ("validSetSize",)

class DrugExValidationStrategyInitSerializer(DrugExValidationStrategySerializer):
    metrics = serializers.PrimaryKeyRelatedField(many=True, queryset=ModelPerformanceMetric.objects.all())

    class Meta:
        model = models.DrugExValidationStrategy
        fields = DrugExValidationStrategySerializer.Meta.fields + ("validSetSize",)

class DrugExTrainingStrategySerializer(TrainingStrategySerializer):

    class Meta:
        model = models.DrugExTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields

class DrugExTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):

    class Meta:
        model = models.DrugExTrainingStrategy
        fields = TrainingStrategyInitSerializer.Meta.fields


class DrugExNetSerializer(ModelSerializer):
    molset = MolSetSerializer(many=False)
    corpus = DrugExCorpusSerializer(many=False)
    trainingStrategy = DrugExTrainingStrategySerializer(many=False)
    validationStrategy = DrugExValidationStrategySerializer(many=False)

    class Meta:
        model = models.DrugExNet
        fields = ModelSerializer.Meta.fields + ("molset", "corpus")
        read_only_fields = ModelSerializer.Meta.read_only_fields + ("corpus",)


class DrugExNetInitSerializer(DrugExNetSerializer):
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=MolSet.objects.all(), required=True)
    trainingStrategy = DrugExTrainingStrategyInitSerializer(many=False)
    validationStrategy = DrugExValidationStrategyInitSerializer(many=False)

    class Meta:
        model = models.DrugExNet
        fields = ("name", "description", "project", "molset", "trainingStrategy", "validationStrategy")

    def create(self, validated_data):
        instance = models.DrugExNet.objects.create(
            name=validated_data['name'],
            description=validated_data['description'],
            project=validated_data['project'],
            molset=validated_data['molset'],
        )

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.DrugExTrainingStrategy.objects.create(
            modelInstance=instance,
            algorithm=strat_data['algorithm'],
            mode=strat_data['mode'],
        )

        self.saveParameters(trainingStrategy, strat_data)

        strat_data = validated_data['validationStrategy']
        validationStrategy = models.DrugExValidationStrategy.objects.create(
            modelInstance = instance,
            validSetSize=strat_data['validSetSize']
        )
        validationStrategy.metrics.set(strat_data['metrics'])
        validationStrategy.save()

        return instance