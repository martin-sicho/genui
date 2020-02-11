"""
serializers

Created by: Martin Sicho
On: 13-01-20, 11:07
"""

from rest_framework import serializers

import modelling.models
from compounds.serializers import MolSetSerializer
from modelling.serializers import TrainingStrategySerializer, BasicValidationStrategyInitSerializer, ModelSerializer, \
    BasicValidationStrategySerializer, TrainingStrategyInitSerializer
from . import models


class DescriptorGroupSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = models.DescriptorGroup
        fields = ('id', 'name',)

class QSARTrainingStrategySerializer(TrainingStrategySerializer):
    descriptors = DescriptorGroupSerializer(many=True)

    class Meta:
        model = models.QSARTrainingStrategy
        fields = TrainingStrategySerializer.Meta.fields + ('descriptors', 'activityThreshold')

class QSARTrainingStrategyInitSerializer(TrainingStrategyInitSerializer):
    descriptors = serializers.PrimaryKeyRelatedField(many=True, queryset=models.DescriptorGroup.objects.all())

    class Meta:
        model = models.QSARTrainingStrategy
        fields = QSARTrainingStrategySerializer.Meta.fields


class QSARModelSerializer(ModelSerializer):
    trainingStrategy = QSARTrainingStrategySerializer(many=False)
    validationStrategy = BasicValidationStrategySerializer(many=False)
    molset = MolSetSerializer(many=False)
    predictions = serializers.PrimaryKeyRelatedField(many=True, queryset=models.ActivitySet.objects.all())
    taskID = serializers.UUIDField(required=False)

    class Meta:
        model = models.QSARModel
        fields = ModelSerializer.Meta.fields + ('molset', 'predictions', 'taskID')
        read_only_fields = ModelSerializer.Meta.read_only_fields + ('predictions', 'taskID')

class QSARModelInitSerializer(QSARModelSerializer):
    trainingStrategy = QSARTrainingStrategyInitSerializer(many=False)
    validationStrategy = BasicValidationStrategyInitSerializer(many=False)
    molset = serializers.PrimaryKeyRelatedField(many=False, queryset=models.MolSet.objects.all())

    class Meta:
        model = models.QSARModel
        fields = ('name', 'description', 'project', 'trainingStrategy', 'validationStrategy', 'molset')
        read_only_fields = None

    def create(self, validated_data, **kwargs):
        instance = super().create(
            validated_data
            , molset=validated_data['molset']
            , **kwargs
        )
        models.ModelActivitySet.objects.create(model=instance, project=instance.project)

        strat_data = validated_data['trainingStrategy']
        trainingStrategy = models.QSARTrainingStrategy(
            modelInstance = instance,
            algorithm = strat_data['algorithm'],
            mode = strat_data['mode'],
            activityThreshold = strat_data['activityThreshold']
        )
        trainingStrategy.save()
        trainingStrategy.descriptors.set(strat_data['descriptors'])
        trainingStrategy.save()

        self.saveParameters(trainingStrategy, strat_data)

        strat_data = validated_data['validationStrategy']
        validationStrategy = modelling.models.BasicValidationStrategy.objects.create(
            modelInstance = instance,
            cvFolds=strat_data['cvFolds'],
            validSetSize=strat_data['validSetSize']
        )
        validationStrategy.metrics.set(strat_data['metrics'])
        validationStrategy.save()

        return instance

